# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

"""
Helper functions for constructing and filtering VariantSequence objects
from reads overlapping a variant locus.
"""

from bisect import bisect_left
from collections import defaultdict

import numpy as np

from .allele_read_helpers import group_unique_sequences
from .assembly import merge_identical_sequences
from .logging import get_logger
from .variant_sequence import VariantSequence
from .read_identity import partition_read_alignments, read_sort_key

logger = get_logger(__name__)


def _group_unique_sequences_and_transcript_paths(
        allele_reads,
        max_prefix_size=None,
        max_suffix_size=None):
    """Group identical cDNA observations without merging transcript paths."""
    sequence_groups = group_unique_sequences(
        allele_reads,
        max_prefix_size=max_prefix_size,
        max_suffix_size=max_suffix_size)
    groups = defaultdict(set)
    for (prefix, alt, suffix), reads in sequence_groups.items():
        for read in reads:
            groups[
                prefix,
                alt,
                suffix,
                read.compatible_transcript_ids,
            ].add(read)
    return groups


def initial_variant_sequences_from_reads(
        variant_reads,
        max_nucleotides_before_variant=None,
        max_nucleotides_after_variant=None):
    """
    Get all unique sequences from reads spanning a variant locus. This will
    include partial sequences due to reads starting in the middle of the
    sequence around around a variant.

    Parameters
    ----------
    variant_reads : list of AlleleRead objects

    max_nucleotides_before_variant : int or None

    max_nucleotides_after_variant : int or None

    Returns
    -------
    list of VariantSequence
    """
    unique_sequence_groups = _group_unique_sequences_and_transcript_paths(
        variant_reads,
        max_prefix_size=max_nucleotides_before_variant,
        max_suffix_size=max_nucleotides_after_variant)

    return [
        VariantSequence(
            prefix=prefix,
            alt=alt,
            suffix=suffix,
            reads=compatible_reads)
        for ((prefix, alt, suffix, _compatible_transcript_ids), reads)
        in unique_sequence_groups.items()
        for compatible_reads in partition_read_alignments(reads)
    ]


class _FlankIndex:
    """Linear-space index of strings compatible with an anchored flank.

    Longer matches occupy one contiguous lexicographic range. Shorter matches
    are exact ancestors of the query. Reverse upstream flanks so that both
    sides use prefix matching. Store ranges, not a quadratic table of matches.
    """

    def __init__(self, groups, upstream):
        self.entries = sorted(
            ((prefix[::-1] if upstream else suffix, group_id)
             for group_id, (prefix, suffix, reads, transcript_ids) in enumerate(groups)))
        self.words = [word for word, _ in self.entries]
        self.exact = {}
        for i, word in enumerate(self.words):
            start = self.exact[word][0] if word in self.exact else i
            self.exact[word] = (start, i + 1)
        self.lengths = sorted({len(word) for word in self.exact})

    def compatible_ranges(self, query):
        start = bisect_left(self.words, query)
        low, high = start, len(self.words)
        # No alphabet/sentinel assumptions, including the empty query.
        while low < high:
            middle = (low + high) // 2
            if self.words[middle].startswith(query):
                low = middle + 1
            else:
                high = middle
        ranges = [(start, low)]
        for length in self.lengths:
            if length >= len(query):
                break
            interval = self.exact.get(query[:length])
            if interval is not None:
                ranges.append(interval)
        return sum(end - start for start, end in ranges), ranges


def _variant_sequences_with_shared_read_support(variant_sequences, read_groups):
    """Assign compatible original reads to each retained candidate.

    Compatibility is checked against the reads themselves, not transitively
    through other candidates: a short core can match two conflicting longer
    branches without either branch supporting the other. Overhangs outside a
    candidate are irrelevant, and coverage still uses each read's real bounds.
    No sequence is extended or invented here. Reads may support more than one
    alternative; downstream aggregation deduplicates them by set union.
    ``read_groups`` is the canonical partition from ``group_unique_sequences``.
    """
    variant_sequences = list(variant_sequences)
    if not variant_sequences:
        return []
    max_prefix = max(len(s.prefix) for s in variant_sequences)
    max_suffix = max(len(s.suffix) for s in variant_sequences)
    normalized_read_groups = defaultdict(set)
    for key, reads in read_groups.items():
        if len(key) == 4:
            normalized_read_groups[key].update(reads)
            continue
        prefix, alt, suffix = key
        for read in reads:
            normalized_read_groups[
                prefix,
                alt,
                suffix,
                getattr(read, "compatible_transcript_ids", None),
            ].add(read)
    bounded_groups = defaultdict(set)
    for (prefix, alt, suffix, transcript_ids), reads in normalized_read_groups.items():
        # Only index the observable context; retain the ORIGINAL read objects
        # and their real bounds for coverage. Distinct overhangs outside every
        # candidate can safely share an index entry, never a read identity.
        key = (
            prefix[-max_prefix:] if max_prefix else "",
            alt,
            suffix[:max_suffix],
            transcript_ids)
        bounded_groups[key].update(reads)
    by_allele = defaultdict(list)
    def group_key(item):
        prefix, alt, suffix, ids = item[0]
        return prefix, alt, suffix, (0, ()) if ids is None else (1, tuple(sorted(ids)))

    for (prefix, alt, suffix, transcript_ids), reads in sorted(bounded_groups.items(), key=group_key):
        by_allele[alt].append((prefix, suffix, reads, transcript_ids))
    indexes = {alt: (_FlankIndex(groups, True), _FlankIndex(groups, False))
               for alt, groups in by_allele.items()}
    result = []
    for sequence in variant_sequences:
        supporting_reads = set()
        constraints = dict(sequence._source_alignments)
        boundaries = [0] * (len(sequence) + 1)
        variant_start, variant_end = sequence.variant_indices()
        if sequence.alt in indexes:
            upstream, downstream = indexes[sequence.alt]
            n_left, left_ranges = upstream.compatible_ranges(sequence.prefix[::-1])
            n_right, right_ranges = downstream.compatible_ranges(sequence.suffix)
            index, ranges = ((upstream, left_ranges) if n_left <= n_right else (downstream, right_ranges))
            groups = by_allele[sequence.alt]
            for start, end in ranges:
                for i in range(start, end):
                    prefix, suffix, group, transcript_ids = groups[index.entries[i][1]]
                    # The selective flank has already matched; check both
                    # here to keep the exact compatibility rule explicit.
                    path_supports_candidate = (
                        sequence.compatible_transcript_ids is None
                        and transcript_ids is None
                    ) or (
                        sequence.compatible_transcript_ids is not None
                        and (
                            transcript_ids is None
                            or sequence.compatible_transcript_ids.issubset(
                                transcript_ids)
                        )
                    )
                    if (path_supports_candidate
                            and (sequence.alt or (prefix and sequence.prefix) or (suffix and sequence.suffix))
                            and (prefix.endswith(sequence.prefix) or sequence.prefix.endswith(prefix))
                            and (suffix.startswith(sequence.suffix) or sequence.suffix.startswith(suffix))):
                        compatible_group = []
                        for read in sorted(group, key=read_sort_key):
                            placements = getattr(read, "source_alignments", ())
                            if all(key not in constraints or constraints[key] == value
                                   for key, value in placements):
                                constraints.update(placements)
                                compatible_group.append(read)
                        supporting_reads.update(compatible_group)
                        # Legacy metadata-free groups retain the fast bounded
                        # coverage calculation. Collected observations need
                        # segment-aware coverage below to avoid recounting
                        # a mate already represented in a collapsed pair.
                        boundaries[max(0, variant_start - len(prefix))] += len(compatible_group)
                        boundaries[min(len(sequence), variant_end + len(suffix))] -= len(compatible_group)
        if supporting_reads == sequence.reads:
            supported = sequence
        else:
            supported = VariantSequence(
                prefix=sequence.prefix,
                alt=sequence.alt,
                suffix=sequence.suffix,
                reads=supporting_reads)
        if supported._coverage_cache is None and not constraints:
            supported._coverage_cache = np.cumsum(boundaries, dtype="int32")[:-1]
        result.append(supported)
    return result


def trim_variant_sequences(variant_sequences, min_variant_sequence_coverage):
    """
    Trim VariantSequences to desired coverage and combine exact duplicates.

    Contained but non-identical candidates remain separate because reference
    compatibility is not known at this stage. Reconcile their compatible read
    support before applying the coverage threshold, including when overlap
    assembly is disabled.

    Parameters
    ----------
    variant_sequences : list of VariantSequence

    min_variant_sequence_coverage : int

    Returns list of VariantSequence
    """
    n_total = len(variant_sequences)
    variant_sequences = merge_identical_sequences(
        sequence for sequence in variant_sequences if len(sequence) > 0)
    reads = {read for sequence in variant_sequences for read in sequence.reads}
    # Compare once per distinct read sequence, rather than once per fragment
    # or per appearance in the retained assembly history. Keep this original
    # pool even when some of its candidates are discarded below.
    read_groups = group_unique_sequences(reads)
    variant_sequences = _variant_sequences_with_shared_read_support(
        variant_sequences, read_groups)
    # when no base of a sequence has sufficient coverage, trim_by_coverage
    # returns an empty sequence which still carries every one of the original
    # reads. Those sequences are useless downstream and their read counts are
    # actively misleading when comparing candidates by support, so drop them
    # here rather than relying on a later validity check to sweep them up.
    trimmed_variant_sequences = [
        trimmed
        for trimmed in (
            variant_sequence.trim_by_coverage(min_variant_sequence_coverage)
            for variant_sequence in variant_sequences
        )
        if len(trimmed) > 0
    ]
    trimmed_variant_sequences = merge_identical_sequences(
        trimmed_variant_sequences)
    # Removing low-coverage overhangs can make previously conflicting reads
    # compatible. Recount against the original pool without extending the
    # trimmed sequences or introducing another round of candidate generation.
    trimmed_variant_sequences = _variant_sequences_with_shared_read_support(
        trimmed_variant_sequences, read_groups)
    n_after_trimming = len(trimmed_variant_sequences)
    logger.info(
        "Kept %d/%d variant sequences after read coverage trimming to >=%dx",
        n_after_trimming,
        n_total,
        min_variant_sequence_coverage)
    return trimmed_variant_sequences
