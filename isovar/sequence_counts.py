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

from collections import Counter, defaultdict, namedtuple

from .common import group_unique_sequences, get_variant_nucleotides


SequenceCounts = namedtuple(
    "SequenceCounts",
    [
        # dictionary mapping unique sequences to the sum of the number
        # of reads fully supporting that sequence and fraction of that
        # sequence supported by partially overlapping reads
        "sequence_weights",
        # dictionary of context sequences mapping to the # of reads
        # they were detected in
        "fully_supporting_read_counts",
        # dictionary of context sequences mapping to the # of reads
        # which partially overlap them in (and agree at all overlapped bases)
        "partially_supporting_read_counts",
        # dictionary mapping context sequences to the sum of
        # overlap weights from all partially overlapping reads
        "partially_supporting_read_weights",
        # variant nucleotide string, since all the dictionary keys
        # above are only (prefix, suffix) pairs which are missing these
        # nucleotides
        "variant_nucleotides",
    ])

def sequence_counts(
        variant_reads,
        context_size=45,
        min_sequence_length=None):
    """
    Returns a dictionary mapping (prefix, variant, suffix) to a pair of
    integers: the first indicating the number of reads which fully support
    the sequence and the second is the number of reads which partially
    support the sequence (weighted by the fraction of nucleotides
    overlapping).
    """

    # Get all unique sequences from reads spanning the
    # variant locus. This will include partial sequences
    # due to reads starting in the middle of the context,
    # we'll use these only to compute partial support for a
    # full length sequence
    unique_sequence_groups = group_unique_sequences(
        variant_reads,
        max_prefix_size=context_size,
        max_suffix_size=context_size)

    variant_seq = get_variant_nucleotides(variant_reads)
    variant_len = len(variant_seq)

    if not min_sequence_length:
        # if not specified, then use the longest sequence available
        # which will be at most twice the context size
        min_sequence_length = variant_len + max(
            len(prefix) + len(suffix)
            for (prefix, suffix) in unique_sequence_groups.keys())

    full_sequences = {
        (prefix, suffix): read_names
        for ((prefix, suffix), read_names)
        in unique_sequence_groups.items()
        if (len(prefix) + len(suffix) + variant_len) >= min_sequence_length
    }

    full_sequence_counts = {
        key: len(read_names)
        for (key, read_names)
        in full_sequences.items()
    }

    # index the largest sequences by prefix and suffix, so that it's
    # easy to determine when one sequence is contained within another
    prefix_to_suffixes = defaultdict(set)
    for (prefix, suffix) in full_sequences.keys():
        prefix_to_suffixes[prefix].add(suffix)

    partially_supporting_read_counts = Counter()
    partially_supporting_read_weights = Counter()
    for (prefix, suffix), reads in unique_sequence_groups.items():
        n_reads = len(reads)
        curr_len = len(prefix) + len(suffix)
        for other_prefix, other_suffixes in prefix_to_suffixes.items():
            if other_prefix.endswith(prefix):
                for other_suffix in other_suffixes:
                    if prefix == other_prefix and suffix == other_suffix:
                        # we've already counter exact matches
                        continue
                    elif other_suffix.startswith(suffix):
                        other_len = len(other_prefix) + len(other_suffix)
                        other_key = (other_prefix, other_suffix)
                        fraction_bases_covered = float(curr_len) / other_len
                        assert fraction_bases_covered < 1.0
                        partially_supporting_read_counts[other_key] += n_reads
                        partially_supporting_read_weights[other_key] += (
                            n_reads * fraction_bases_covered)

    combined_weights = {
        key: full_count + partially_supporting_read_weights[key]
        for (key, full_count) in full_sequence_counts.items()
    }
    return SequenceCounts(
        sequence_weights=combined_weights,
        fully_supporting_read_counts=full_sequence_counts,
        partially_supporting_read_counts=partially_supporting_read_counts,
        partially_supporting_read_weights=partially_supporting_read_weights,
        variant_nucleotides=variant_seq)
