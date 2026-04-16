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
Greedy overlap assembly of VariantSequence objects.

This module assembles short cDNA sequences (centered on a variant locus) into
longer contiguous sequences by iteratively merging pairs with the largest
overlaps.

Assumptions / limitations:

- **Exact sequence match**: overlap detection requires exact prefix/suffix
  containment (str.endswith / str.startswith). A single sequencing error in
  the flanking sequence prevents a merge. This is suitable for Illumina
  short reads (~0.1% error rate) but will produce fragmented assemblies
  with long-read technologies (PacBio/ONT, ~5-10% indel error rate).

- **Shared alt allele**: all input VariantSequence objects must carry the
  same alt allele string. Reads with different alleles at the variant locus
  should be separated before assembly.
"""

from collections import defaultdict

from .default_parameters import MIN_VARIANT_SEQUENCE_ASSEMBLY_OVERLAP_SIZE
from .logging import get_logger

logger = get_logger(__name__)


def greedy_merge_helper(
        variant_sequences,
        min_overlap_size=MIN_VARIANT_SEQUENCE_ASSEMBLY_OVERLAP_SIZE):
    """
    Score all valid pairwise merges by overlap size and greedily accept
    the best non-conflicting merges (each sequence participates in at most
    one merge per round). This makes the output deterministic regardless
    of input order.

    Returns a list of merged VariantSequence objects, and True if any
    were successfully merged.
    """
    candidates = []
    for i in range(len(variant_sequences)):
        for j in range(i + 1, len(variant_sequences)):
            combined = variant_sequences[i].combine(
                variant_sequences[j], min_overlap_size=min_overlap_size)
            if combined is not None:
                overlap = (
                    len(variant_sequences[i])
                    + len(variant_sequences[j])
                    - len(combined))
                candidates.append((overlap, len(combined.reads), i, j, combined))

    if not candidates:
        return list(variant_sequences), False

    # largest overlap first, then most reads, then indices for determinism
    candidates.sort(key=lambda c: (-c[0], -c[1], c[2], c[3]))

    used = set()
    merged_variant_sequences = {}
    for overlap, n_reads, i, j, combined in candidates:
        if i in used or j in used:
            continue
        used.add(i)
        used.add(j)
        if combined.sequence in merged_variant_sequences:
            existing = merged_variant_sequences[combined.sequence]
            combined = combined.add_reads(existing.reads)
        merged_variant_sequences[combined.sequence] = combined

    unmerged = [
        variant_sequences[k]
        for k in range(len(variant_sequences))
        if k not in used
    ]
    result = list(merged_variant_sequences.values()) + unmerged
    return result, True


def greedy_merge(
        variant_sequences,
        min_overlap_size=MIN_VARIANT_SEQUENCE_ASSEMBLY_OVERLAP_SIZE):
    """
    Greedily merge overlapping sequences into longer sequences.

    Accepts a collection of VariantSequence objects and returns another
    collection of elongated variant sequences. The reads field of the
    returned VariantSequence object will contain reads which
    only partially overlap the full sequence.
    """
    merged_any = True
    while merged_any:
        variant_sequences, merged_any = greedy_merge_helper(
            variant_sequences,
            min_overlap_size=min_overlap_size)
    return variant_sequences


def collapse_substrings(variant_sequences):
    """
    Combine shorter sequences which are fully contained in longer sequences.

    Parameters
    ----------
    variant_sequences : list
       List of VariantSequence objects

    Returns a (potentially shorter) list without any contained subsequences.
    """
    if len(variant_sequences) <= 1:
        # if we don't have at least two VariantSequences then just
        # return your input
        return variant_sequences

    # dictionary mapping VariantSequence objects to lists of reads
    # they absorb from substring VariantSequences
    extra_reads_from_substrings = defaultdict(set)
    result_list = []
    # sort by longest to shortest total length
    for short_variant_sequence in sorted(
            variant_sequences,
            key=lambda seq: -len(seq)):
        found_superstring = False
        for long_variant_sequence in result_list:
            if long_variant_sequence.contains(short_variant_sequence):
                extra_reads_from_substrings[long_variant_sequence].update(
                    short_variant_sequence.reads)
                found_superstring = True
                break
        if not found_superstring:
            result_list.append(short_variant_sequence)
    # add to each VariantSequence the reads it absorbed from dropped substrings
    # and then return
    return [
        variant_sequence.add_reads(
            extra_reads_from_substrings[variant_sequence])
        for variant_sequence in result_list
    ]


DEFAULT_MAX_ASSEMBLY_SEQUENCES = 1000


def iterative_overlap_assembly(
        variant_sequences,
        min_overlap_size=MIN_VARIANT_SEQUENCE_ASSEMBLY_OVERLAP_SIZE,
        max_assembly_sequences=DEFAULT_MAX_ASSEMBLY_SEQUENCES):
    """
    Assembles longer sequences from reads centered on a variant by
    merging all pairs of overlapping sequences and collapsing
    shorter sequences onto every longer sequence which contains them.

    Parameters
    ----------
    variant_sequences : list of VariantSequence

    min_overlap_size : int

    max_assembly_sequences : int or None
        If the number of input sequences exceeds this threshold after
        substring collapse, skip the O(n^2)-per-round greedy merge and
        return sequences sorted by read support. Set to None to disable.

    Returns a list of variant sequences, sorted by decreasing read support.
    """
    if len(variant_sequences) <= 1:
        return variant_sequences

    n_before_collapse = len(variant_sequences)
    variant_sequences = collapse_substrings(variant_sequences)
    n_after_collapse = len(variant_sequences)
    logger.info(
        "Collapsed %d -> %d sequences",
        n_before_collapse,
        n_after_collapse)

    if (max_assembly_sequences is not None
            and n_after_collapse > max_assembly_sequences):
        logger.warning(
            "Too many variant sequences (%d > %d) for greedy assembly; "
            "skipping merge and returning sequences sorted by read support",
            n_after_collapse,
            max_assembly_sequences)
    else:
        variant_sequences = greedy_merge(variant_sequences, min_overlap_size)

    return list(sorted(
        variant_sequences,
        key=lambda seq: -len(seq.reads)))
