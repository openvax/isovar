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

def one_step_assembly(pair_counts):
    # matrix from seq1->seq2->count
    # where the count indicates the number of reads that support
    # this overlap

    greedy_merges = {}

    def sort_by_decreasing_prefix_length(x):
        return -len(x[0][0])

    def sort_by_decreasing_suffix_length(x):
        return -len(x[0][1])

    def sort_by_increasing_total_length(x):
        return (len(x[0][0]) + len(x[0][1]))

    def sort_by_decreasing_total_length(x):
        return -(len(x[0][0]) + len(x[0][1]))

    # first, greedily merge long prefixes with long suffixes

    for (p1, s1), c1 in sorted(
            pair_counts.items(),
            key=sort_by_decreasing_prefix_length):
        print("p1,s1,c1", p1, s1, c1)
        for (p2, s2), c2 in sorted(
                pair_counts.items(),
                key=sort_by_decreasing_suffix_length):
            if len(p2) > len(p1):
                # only want to merge cases where first sequence is a
                # prefix of the second sequence
                break
            elif p1 == p2 and s1 == s2:
                break
            elif p1.endswith(p2) and s2.startswith(s1):
                # is the candidate sequence is a prefix of the accepted?
                # Example:
                # p1 s1 = XXXXXXXX ZZZZZZ
                # p2 s2 =       XX ZZZZZZZZZ
                # ...
                # then combine them into a longer sequence
                longer_sequence_pair = (p1, s2)
                n_combinations = c1 + c2
                if longer_sequence_pair in greedy_merges:
                    greedy_merges[longer_sequence_pair] += n_combinations
                else:
                    greedy_merges[longer_sequence_pair] = n_combinations
    return greedy_merges
    """
    # now, clean up the set of merged reads by eliminating shorter sequences
    # that are contained in the longer ones

    simplified_greedy_merges = {}
    for (p1, s1), c1 in sorted(
        greedy_merges.items(),
        key=)
    return assembly_counts
    """
def recursive_assembly(prefix_suffix_pairs, n=1):
    assembly_counts = unique_counts(prefix_suffix_pairs)
    for _ in range(n):
        assembly_counts = one_step_assembly(assembly_counts)
    return assembly_counts

def assemble_transcript_fragments(partitioned_read_sequences):
    variant_seq, prefix_suffix_pairs = drop_variant_from_partitioned_sequences(
        partitioned_read_sequences)

    assembly_counts = recursive_assembly(prefix_suffix_pairs)
    return [
        ((prefix, variant_seq, suffix), count)
        for ((prefix, suffix), count)
        in sorted(assembly_counts.items(), key=lambda x: (x[1], len(x[0][0]) + len(x[0][1])))
    ]
