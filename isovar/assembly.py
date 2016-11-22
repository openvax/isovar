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
import logging

from .default_parameters import MIN_VARIANT_SEQUENCE_ASSEMBLY_OVERLAP_SIZE

logger = logging.getLogger(__name__)

def sort_by_decreasing_prefix_length(seq):
    """
    Key function for sorting from longest to shortest prefix length.

    Parameters
    ----------
    seq : VariantSequence
    """
    return -len(seq.prefix)

def sort_by_decreasing_suffix_length(seq):
    """
    Key function for sorting from longest to shortest suffix length.

    Parameters
    ----------
    seq : VariantSequence
    """
    return -len(seq.suffix)

def sort_by_increasing_total_length(seq):
    """
    Key function for sorting from shortest to longest total length.

    Parameters
    ----------
    seq : VariantSequence
    """
    return len(seq.sequence)

def sort_by_decreasing_total_length(seq):
    """
    Key function for sorting from longest to shortest total length.

    Parameters
    ----------
    seq : VariantSequence
    """
    return -len(seq.sequence)

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

    # dictionary mapping distinct sequences to
    # combined VariantSequence objects
    merged_variant_sequences = {}

    # The basic idea here is that we're comparing variant sequences with
    # the longest prefix against those with the longest suffix and trying
    # to combine them into even longer sequences.
    for variant_sequence1 in sorted(
            variant_sequences,
            key=sort_by_decreasing_prefix_length):
        # rather than generating candidate VariantSequence objects for all
        # pairs of original sequences we check each sequence against the
        # new sequences generated thus far. If any of them contain the
        # sequence we're currently considering then just merge it into
        # those sequences instead.
        contained_in_existing_supersequence = False
        for key, candidate_supersequence in list(merged_variant_sequences.items()):
            if candidate_supersequence.contains(variant_sequence1):
                contained_in_existing_supersequence = True
                merged_variant_sequences[key] = candidate_supersequence.combine(variant_sequence1)

        if contained_in_existing_supersequence:
            continue

        # If this variant sequence wasn't contained in any existing
        # merged supersequence, try merging it with all the other inputs.
        # Either we'll find something to merge it with or we'll add the
        # sequence on its own.
        found_merge_partner = False

        for variant_sequence2 in sorted(
                variant_sequences,
                key=sort_by_decreasing_suffix_length):
            if not variant_sequence1.left_overlaps(
                    variant_sequence2,
                    min_overlap_size=min_overlap_size):
                continue
            found_merge_partner = True
            combined = variant_sequence1.combine(variant_sequence2)
            if combined.sequence not in merged_variant_sequences:
                merged_variant_sequences[combined.sequence] = combined
            else:
                # it's possible to get the same merged sequence from distinct
                # input sequences
                # For example
                #   abcXYZddd + cXYZdddd = abcXYZdddd
                #   abcXYZd + bcXYZdddd =  abcXYZdddd
                # In this case we make a VariantSequence record with the
                # reads from both original sequences.
                existing_record_with_same_sequence = merged_variant_sequences[
                    combined.sequence]
                combined_with_more_reads = existing_record_with_same_sequence.combine(
                    combined)
                merged_variant_sequences[combined.sequence] = combined_with_more_reads
        if not found_merge_partner:
            # if we weren't able to merge this sequence with any other
            # sequences then add it by itself
            merged_variant_sequences[variant_sequence1.sequence] = variant_sequence1
    return list(merged_variant_sequences.values())

def collapse_substrings(variant_sequences):
    """
    Combine shorter sequences which are fully contained in longer sequences.

    Parameters
    ----------
    variant_sequences : list
       List of VariantSequence objects

    Returns a (potentially shorter) list without any contained subsequences.
    """
    # dictionary mapping VariantSequence objects to lists of reads
    # they absorb from substring VariantSequences
    extra_reads_from_substrings = defaultdict(set)
    result_list = []
    for short_variant_sequence in sorted(
            variant_sequences,
            key=sort_by_decreasing_total_length):
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
        variant_sequence.add_reads(extra_reads_from_substrings[variant_sequence])
        for variant_sequence in result_list
    ]

def sort_by_decreasing_read_count_and_sequence_lenth(variant_sequence):
    """
    Sort variant sequences by number of supporting reads and length of
    assembled sequence.
    """
    return -len(variant_sequence.reads), -len(variant_sequence.sequence)

def iterative_overlap_assembly(
        variant_sequences,
        min_overlap_size=MIN_VARIANT_SEQUENCE_ASSEMBLY_OVERLAP_SIZE,
        n_merge_iters=2):
    """
    Assembles longer sequences from reads centered on a variant by
    between merging all pairs of overlapping sequences and collapsing
    shorter sequences onto every longer sequence which contains them.
    """
    # reduce the number of inputs to the merge algorithm by first collapsing
    # shorter sequences onto the longer sequences which contain them
    n_before_collapse = len(variant_sequences)

    variant_sequences = collapse_substrings(variant_sequences)
    n_after_collapse = len(variant_sequences)
    logger.info(
        "Collapsed %d -> %d sequences",
        n_before_collapse,
        n_after_collapse)

    for i in range(n_merge_iters):
        n_before_merge = len(variant_sequences)
        variant_sequences_after_merge = greedy_merge(
            variant_sequences,
            min_overlap_size=min_overlap_size)
        n_after_merge = len(variant_sequences_after_merge)
        logger.info(
            "Assembly iter %d/%d: merged %d VariantSequences into %d",
            i + 1,
            n_merge_iters,
            n_before_merge,
            n_after_merge)
        # did the set of variant sequences change? if not, then we're done
        if set(vs.sequence for vs in variant_sequences) == set(
                vs.sequence for vs in variant_sequences_after_merge):
            logger.info("Converged on iter #%d with %d sequences" % (
                i + 1, n_after_collapse))
            break

        if n_after_merge == 0:
            # if the greedy merge procedure fails for all pairs of candidate
            # sequences then we'll get an empty set of new longer sequences,
            # in which case we should just stop with the results of the last
            # iteration
            return variant_sequences
        elif n_after_merge == 1:
            # once we have only one sequence then there's no point trying
            # to further combine sequences
            return variant_sequences_after_merge
        variant_sequences = variant_sequences_after_merge
    # final cleanup step, merge any VariantSequences which contain each other
    return collapse_substrings(variant_sequences)
