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

from .common import group_unique_sequences, get_variant_nucleotides


def sort_by_decreasing_prefix_length(x):
    return -len(x[0][0])

def sort_by_decreasing_suffix_length(x):
    return -len(x[0][1])

def sort_by_increasing_total_length(x):
    return (len(x[0][0]) + len(x[0][1]))

def sort_by_decreasing_total_length(x):
    return -(len(x[0][0]) + len(x[0][1]))

def greedy_merge(prefix_suffix_pair_groups, min_overlap_size=30):
    """
    Greedily merge overlapping sequences into longer sequences.

    Returns a dictionary mapping (prefix, suffix) sequence pairs to
    counts of how many reads support the new longer sequences.
    """
    merged = {}
    for (p1, s1), names1 in sorted(
            prefix_suffix_pair_groups.items(),
            key=sort_by_decreasing_prefix_length):
        len1 = len(p1) + len(s1)
        for (p2, s2), names2 in sorted(
                prefix_suffix_pair_groups.items(),
                key=sort_by_decreasing_suffix_length):
            if len(p2) > len(p1):
                continue
            elif len(s2) < len(s1):
                continue
            possible_overlap = len(p2) + len(s1)

            if possible_overlap < min_overlap_size:
                continue

            len2 = len(p2) + len(s2)
            new_length = len(p1) + len(s2)
            # if we're not extending the sequence length, then why bother?
            if new_length <= len1 or new_length <= len2:
                continue
            elif p1.endswith(p2) and s2.startswith(s1):
                # is the candidate sequence is a prefix of the accepted?
                # Example:
                # p1 s1 = XXXXXXXX ZZZZZZ
                # p2 s2 =       XX ZZZZZZZZZ
                # ...
                # then combine them into a longer sequence
                longer_pair = (p1, s2)
                combined_read_names = names1.union(names2)
                if longer_pair in merged:
                    merged[longer_pair] = merged[longer_pair].union(
                        combined_read_names)
                else:
                    merged[longer_pair] = combined_read_names
    return merged

def collapse_substrings(assembly_groups):
    """
    Combine shorter sequences which are fully contained in longer sequences
    """
    assert len(assembly_groups) > 0
    sorted_pairs = list(sorted(
        assembly_groups.items(),
        key=sort_by_decreasing_total_length))
    result_list = []
    for ((prefix_short, suffix_short), names_short) in sorted_pairs:
        found_superstring = False
        for (prefix_long, suffix_long), names_long in result_list:
            if prefix_long.endswith(prefix_short) and suffix_long.startswith(
                    suffix_short):
                found_superstring = True
                for name in names_short:
                    names_long.add(name)
        if not found_superstring:
            result_list.append(
                ((prefix_short, suffix_short), names_short.copy()))
    return {
        (p, s): names
        for ((p, s), names) in result_list
    }

def recursive_assembly(variant_reads, min_overlap_size=30, n_merge_iters=2):
    """
    Assembles longer sequences from reads centered on a variant by alternating
    between merging all pairs of overlapping sequences and collapsing
    shorter sequences onto every longer sequence which contains them.
    """
    assert len(variant_reads) > 0
    assembly_groups = group_unique_sequences(variant_reads)
    for i in range(n_merge_iters):
        previous_assembly_groups = {k: v for (k, v) in assembly_groups.items()}

        assembly_groups = greedy_merge(
            assembly_groups,
            min_overlap_size=min_overlap_size)

        if len(assembly_groups) == 0:
            return previous_assembly_groups
        assembly_groups = collapse_substrings(assembly_groups)
        if len(assembly_groups) == 1:
            return assembly_groups
    return assembly_groups

def final_sort_key(x):
    supporting_read_names = x[1]
    n_supporting_reads = len(supporting_read_names)
    prefix = x[0][0]
    suffix = x[0][1]
    combined_length = len(prefix) + len(suffix)
    return (n_supporting_reads, combined_length)

def assemble_transcript_fragments(
        variant_reads,
        min_overlap_size=30,
        n_merge_iters=2):
    assert len(variant_reads) > 0
    variant_seq = get_variant_nucleotides(variant_reads)
    assembly_counts = recursive_assembly(
        variant_reads,
        min_overlap_size=min_overlap_size,
        n_merge_iters=n_merge_iters)
    return [
        ((prefix, variant_seq, suffix), len(supporting_read_names))
        for ((prefix, suffix), supporting_read_names)
        in sorted(assembly_counts.items(), key=final_sort_key)
    ]
