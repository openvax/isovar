# Copyright (c) 2016. Mount Sinai School of Medicine
#
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

from __future__ import print_function, division, absolute_import

from collections import defaultdict

from six.moves import range

from .default_parameters import MIN_VARIANT_SEQUENCE_ASSEMBLY_OVERLAP_SIZE
from .logging import get_logger

logger = get_logger(__name__)


def greedy_merge_helper(
        variant_sequences,
        min_overlap_size=MIN_VARIANT_SEQUENCE_ASSEMBLY_OVERLAP_SIZE):
    """
    Returns a list of merged VariantSequence objects, and True if any
    were successfully merged.
    """
    merged_variant_sequences = {}
    merged_any = False

    # here we'll keep track of sequences that haven't been merged yet, and add them in at the end
    unmerged_variant_sequences = set(variant_sequences)
    for i in range(len(variant_sequences)):
        sequence1 = variant_sequences[i]
        # it works to loop over the triangle (i+1 onwards) because combine() tries flipping the
        # arguments if sequence1 is on the right of sequence2
        for j in range(i + 1, len(variant_sequences)):
            sequence2 = variant_sequences[j]
            combined = sequence1.combine(sequence2)
            if combined is None:
                continue
            if combined.sequence in merged_variant_sequences:
                existing = merged_variant_sequences[combined.sequence]
                # the existing VariantSequence and the newly merged
                # VariantSequence should differ only in which reads support them
                combined = combined.add_reads(existing.reads)
            merged_variant_sequences[combined.sequence] = combined
            unmerged_variant_sequences.discard(sequence1)
            unmerged_variant_sequences.discard(sequence2)
            merged_any = True
    result = list(merged_variant_sequences.values()) + list(unmerged_variant_sequences)
    return result, merged_any

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
            found_superstring = long_variant_sequence.contains(short_variant_sequence)
            if found_superstring:
                extra_reads_from_substrings[long_variant_sequence].update(
                    short_variant_sequence.reads)
        if not found_superstring:
            result_list.append(short_variant_sequence)
    # add to each VariantSequence the reads it absorbed from dropped substrings
    # and then return
    return [
        variant_sequence.add_reads(
            extra_reads_from_substrings[variant_sequence])
        for variant_sequence in result_list
    ]

def iterative_overlap_assembly(
        variant_sequences,
        min_overlap_size=MIN_VARIANT_SEQUENCE_ASSEMBLY_OVERLAP_SIZE):
    """
    Assembles longer sequences from reads centered on a variant by
    between merging all pairs of overlapping sequences and collapsing
    shorter sequences onto every longer sequence which contains them.

    Returns a list of variant sequences, sorted by decreasing read support.
    """
    if len(variant_sequences) <= 1:
        # if we don't have at least two sequences to start with then
        # skip the whole mess below
        return variant_sequences

    # reduce the number of inputs to the merge algorithm by first collapsing
    # shorter sequences onto the longer sequences which contain them
    n_before_collapse = len(variant_sequences)
    variant_sequences = collapse_substrings(variant_sequences)
    n_after_collapse = len(variant_sequences)
    logger.info(
        "Collapsed %d -> %d sequences",
        n_before_collapse,
        n_after_collapse)

    merged_variant_sequences = greedy_merge(variant_sequences, min_overlap_size)
    return list(sorted(
        merged_variant_sequences,
        key=lambda seq: -len(seq.reads)))
