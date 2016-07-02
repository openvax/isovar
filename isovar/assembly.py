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

from .variant_sequence import (
    supporting_reads_to_variant_sequences,
    VariantSequence
)

MIN_OVERLAP_SIZE = 30

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

def combine_variant_sequences(
        variant_sequence1,
        variant_sequence2,
        min_overlap_size=MIN_OVERLAP_SIZE):
    """
    Parameters
    ----------
    variant_sequence1 : VariantSequence

    variant_sequence2 : VariantSequence

    Returns either an elongated VariantSequence or None
    """
    if variant_sequence1.alt != variant_sequence2.alt:
            # allele must match!
            return
    if len(variant_sequence2.prefix) > len(variant_sequence1.prefix):
        # only consider strings that overlap like:
        #   variant_sequence1: ppppAssss
        #   variant_sequence2:   ppAsssssss
        # which excludes cases where variant_sequence2 has a longer
        # prefix
        return
    elif len(variant_sequence2.suffix) < len(variant_sequence1.suffix):
        # similarly, we throw cases where variant_sequence2 is shorter
        # after the alt nucleotides than variant_sequence1
        return

    possible_overlap = (
        len(variant_sequence2.prefix) +
        len(variant_sequence1.suffix)
    )

    if possible_overlap < min_overlap_size:
        return

    # compare lengths of the two old sequences and the candidate
    # sequence we're considering constructing
    len1 = len(variant_sequence1.prefix) + len(variant_sequence1.suffix)
    len2 = len(variant_sequence2.prefix) + len(variant_sequence2.suffix)
    new_length = len(variant_sequence1.prefix) + len(variant_sequence2.suffix)

    if new_length <= len1 or new_length <= len2:
        # if we're not extending the sequence length, then why bother?
        return

    prefix_ok = variant_sequence1.prefix.endswith(variant_sequence2.prefix)
    suffix_ok = variant_sequence2.suffix.startswith(variant_sequence1.suffix)
    if prefix_ok and suffix_ok:
        # is the candidate sequence is a prefix of the accepted?
        # Example:
        # p1 a1 s1 = XXXXXXXX Y ZZZZZZ
        # p2 a2 s2 =       XX Y ZZZZZZZZZ
        # ...
        # then combine them into a longer sequence
        return VariantSequence(
            prefix=variant_sequence1.prefix,
            alt=variant_sequence1.alt,
            suffix=variant_sequence2.suffix,
            reads=set(
                variant_sequence1.reads).union(set(variant_sequence2.reads)))

def greedy_merge(variant_sequences, min_overlap_size=MIN_OVERLAP_SIZE):
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
        for variant_sequence2 in sorted(
                variant_sequences,
                key=sort_by_decreasing_suffix_length):
            combined = combine_variant_sequences(
                variant_sequence1,
                variant_sequence2,
                min_overlap_size=min_overlap_size)
            if combined is None:
                continue

            key = combined.sequence

            if key in merged_variant_sequences:
                # it's possible to get the same merged sequence from distinct
                # input sequences
                # For example
                #   abcXYZddd + cXYZdddd = abcXYZdddd
                #   abcXYZd + bcXYZdddd =  abcXYZdddd
                # In this case we make a VariantSequence record with the
                # reads from both original sequences.
                # TODO: In the future I'd like to track how much of each total
                # sequence is supported by any one read.
                existing_record_with_same_sequence = merged_variant_sequences[key]
                union_of_reads = set(combined.reads).union(
                    existing_record_with_same_sequence.reads)
                combined_with_more_reads = VariantSequence(
                    prefix=combined.prefix,
                    alt=combined.alt,
                    suffix=combined.suffix,
                    reads=union_of_reads)
                merged_variant_sequences[key] = combined_with_more_reads
            else:
                merged_variant_sequences[key] = combined
    return merged_variant_sequences

def collapse_substrings(variant_sequences_dict):
    """
    Combine shorter sequences which are fully contained in longer sequences.

    Parameters
    ----------
    variant_sequences_dict : dict
        Dictionary mapping each distinct string sequence
        (prefix + alt + suffix) to a VariantSequence object containing the
        supporting reads for that sequence

    Returns a similar dictionary but with substrings subsumed by longer sequences.
    """
    sorted_pairs = list(sorted(
        assembly_groups.items(),
        key=sort_by_decreasing_total_length))
    result_list = []
    for short_variant_sequence in sorted(
            variant_sequences,
            key=sort_by_decreasing_total_length):
        found_superstring = False
        for long_variant_sequence in result_list:
            if short_variant_sequence.alt != long_variant_sequence.alt:
                # if the alleles disagree then can't combine
                continue
            if prefix_long.endswith(prefix_short) and suffix_long.startswith(
                    suffix_short):
                found_superstring = True
                for name in names_short:
                    names_long.add(name)
        if not found_superstring:
            result_list.append(
                ((prefix_short, allele, suffix_short), names_short.copy()))
    return dict(result_list)


def sort_by_decreasing_read_count_and_sequence_lenth(variant_sequence):
    """
    Sort variant sequences by number of supporting reads and length of
    assembled sequence.
    """
    return -len(variant_sequence.reads), -len(variant_sequence.sequence)

def iterative_assembly(
        variant_sequences,
        min_overlap_size=30,
        n_merge_iters=2):
    """
    Assembles longer sequences from reads centered on a variant by alternating
    between merging all pairs of overlapping sequences and collapsing
    shorter sequences onto every longer sequence which contains them.
    """
    for i in range(n_merge_iters):
        previous_sequences = variant_sequences.copy()
        variant_sequences = greedy_merge(
            variant_sequences,
            min_overlap_size=min_overlap_size)

        if len(variant_sequences) == 0:
            # if the greedy merge procedure fails for all pairs of candidate
            # sequences then we'll get an empty set of new longer sequences,
            # in which case we should just stop with the results of the last
            # iteration
            return previous_sequences

        variant_sequences = collapse_substrings(variant_sequences)

        if len(variant_sequences) == 1:
            # once we have only one sequence then there's no point trying
            # to further combine sequences
            break
    return variant_sequences

def assemble_reads_into_variant_sequences(
        variant_reads,
        min_overlap_size=30,
        n_merge_iters=2,
        initial_sequence_length=None):
    """
    Turn a collection of AlleleRead objects into a collection of VariantSequence
    objects, which are returned in order of (# supported reads, sequence length)
    """
    if initial_sequence_length is None:
        # if no length given for initial 'consensus' sequences then
        # just use the shortest read length
        initial_sequence_length = min(
            len(r.prefix + r.allele + r.suffix)
            for r in variant_reads)

    initial_variant_sequences = supporting_reads_to_variant_sequences(
        reads=variant_reads,
        preferred_sequence_length=initial_sequence_length,
        min_reads_supporting_cdna_sequence=1)
    return iterative_assembly(
        initial_variant_sequences,
        min_overlap_size=min_overlap_size,
        n_merge_iters=n_merge_iters)
