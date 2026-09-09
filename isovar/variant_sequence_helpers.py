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

logger = get_logger(__name__)


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
    unique_sequence_groups = group_unique_sequences(
        variant_reads,
        max_prefix_size=max_nucleotides_before_variant,
        max_suffix_size=max_nucleotides_after_variant)

    return [
        VariantSequence(
            prefix=prefix,
            alt=alt,
            suffix=suffix,
            reads=reads)
        for ((prefix, alt, suffix), reads)
        in unique_sequence_groups.items()
    ]


def filter_variant_sequences_by_read_support(
        variant_sequences,
        min_variant_sequence_coverage):
    """
    Filter VariantSequences to only keep those with at least the desired
    level of coverage.

    Parameters
    ----------
    variant_sequences : list of VariantSequence

    min_variant_sequence_coverage : int
        Minimum number of reads which must cover each
        base of a VariantSequence

    Returns
    -------
    list of VariantSequence
    """
    n_total = len(variant_sequences)
    variant_sequences = [
        s
        for s in variant_sequences
        if s.min_coverage() >= min_variant_sequence_coverage
    ]
    n_dropped = n_total - len(variant_sequences)
    if n_dropped > 0:
        logger.info(
            "Dropped %d/%d variant sequences less than %d supporting reads",
            n_dropped,
            n_total,
            min_variant_sequence_coverage)
    return variant_sequences


def filter_variant_sequences_by_length(
        variant_sequences,
        preferred_sequence_length):
    """
    Retain nonempty variant sequences for downstream translation.

    This function historically dropped every sequence shorter than the
    longest observed sequence (capped at ``preferred_sequence_length``). That
    decision cannot be made safely before translation: a longer sequence may
    fail to match every reference context, raw read counts are only the second
    part of the downstream support ordering, and multiple cDNA sequences may
    be grouped into one protein sequence whose combined support beats every
    individual candidate.

    Keep the legacy function name and ``preferred_sequence_length`` argument
    for API compatibility, but defer all candidate ranking until after
    reference matching and translation. Remove zero-length sequences here as
    a defense in depth; ``trim_variant_sequences`` normally removes them where
    they are produced.

    Parameters
    ----------
    variant_sequences : list of VariantSequence

    preferred_sequence_length : int
        Retained for API compatibility. Initial sequence construction already
        uses this value to bound assembled sequence length.

    Returns
    -------
    list of VariantSequence
    """
    n_total = len(variant_sequences)
    nonempty_variant_sequences = [
        s for s in variant_sequences
        if len(s) > 0
    ]
    n_dropped = n_total - len(nonempty_variant_sequences)
    if n_dropped > 0:
        logger.info(
            "Dropped %d/%d empty variant sequences before translation",
            n_dropped,
            n_total)
    return nonempty_variant_sequences


class _FlankIndex:
    """Linear-space index of strings compatible with an anchored flank.

    Longer matches occupy one contiguous lexicographic range. Shorter matches
    are exact ancestors of the query. Reverse upstream flanks so that both
    sides use prefix matching. Store ranges, not a quadratic table of matches.
    """

    def __init__(self, groups, upstream):
        self.entries = sorted(
            ((prefix[::-1] if upstream else suffix, group_id)
             for group_id, (prefix, suffix, reads) in enumerate(groups)))
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
    bounded_groups = defaultdict(set)
    for (prefix, alt, suffix), reads in read_groups.items():
        # Only index the observable context; retain the ORIGINAL read objects
        # and their real bounds for coverage. Distinct overhangs outside every
        # candidate can safely share an index entry, never a read identity.
        key = (prefix[-max_prefix:] if max_prefix else "", alt, suffix[:max_suffix])
        bounded_groups[key].update(reads)
    by_allele = defaultdict(list)
    for (prefix, alt, suffix), reads in bounded_groups.items():
        by_allele[alt].append((prefix, suffix, reads))
    indexes = {alt: (_FlankIndex(groups, True), _FlankIndex(groups, False))
               for alt, groups in by_allele.items()}
    result = []
    for sequence in variant_sequences:
        supporting_reads = set()
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
                    prefix, suffix, group = groups[index.entries[i][1]]
                    # The selective flank has already matched; check both
                    # here to keep the exact compatibility rule explicit.
                    if ((sequence.alt or (prefix and sequence.prefix) or (suffix and sequence.suffix))
                            and (prefix.endswith(sequence.prefix) or sequence.prefix.endswith(prefix))
                            and (suffix.startswith(sequence.suffix) or sequence.suffix.startswith(suffix))):
                        supporting_reads.update(group)
                        # Canonical read groups partition the original read
                        # objects. Their bounded flanks have identical overlap
                        # bounds within this candidate, so count the group
                        # once instead of rescanning every read for coverage.
                        boundaries[max(0, variant_start - len(prefix))] += len(group)
                        boundaries[min(len(sequence), variant_end + len(suffix))] -= len(group)
        if supporting_reads == sequence.reads:
            supported = sequence
        else:
            supported = VariantSequence(
                prefix=sequence.prefix,
                alt=sequence.alt,
                suffix=sequence.suffix,
                reads=supporting_reads)
        if supported._coverage_cache is None:
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


def filter_variant_sequences(
        variant_sequences,
        preferred_sequence_length,
        min_variant_sequence_coverage):
    """
    Trim variant sequences to the requested coverage and discard any which
    contain no usable sequence afterward.

    Candidate ranking is deliberately deferred until after reference matching
    and translation, when transcript mismatches and aggregate protein support
    are available.

    Parameters
    ----------
    variant_sequences : list of VariantSequence

    preferred_sequence_length : int

    min_variant_sequence_coverage : int

    Returns list of VariantSequence
    """
    variant_sequences = trim_variant_sequences(
        variant_sequences, min_variant_sequence_coverage)

    return filter_variant_sequences_by_length(
        variant_sequences=variant_sequences,
        preferred_sequence_length=preferred_sequence_length)
