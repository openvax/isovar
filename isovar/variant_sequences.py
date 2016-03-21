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

from collections import Counter, defaultdict, namedtuple, OrderedDict

import numpy as np
import pandas as pd

from .common import group_unique_sequences, get_variant_nucleotides
from .variant_reads import variant_reads_generator
from .logging import create_logger

logger = create_logger(__name__)

DEFAULT_SEQUENCE_LENGTH = 105
DEFAULT_CONTEXT_SIZE = DEFAULT_SEQUENCE_LENGTH // 2
DEFAULT_MIN_READS = 1

variant_sequences_object_fields = [
    # dictionary mapping unique sequences to the sum of the number
    # of reads fully supporting that sequence and fraction of that
    # sequence supported by partially overlapping reads
    "combined_sequence_weights",
    # dictionary of context sequences mapping to the # of reads
    # they were detected in
    "full_read_counts",
    # names of reads which fully overlap a sequence
    "full_read_names",
    # dictionary of context sequences mapping to the # of reads
    # which partially overlap them in (and agree at all overlapped bases)
    "partial_read_counts",
    # dictionary from context sequence to the names of partially
    # overlapping reads
    "partial_read_names",
    # dictionary mapping context sequences to the sum of
    # overlap weights from all partially overlapping reads
    "partial_read_weights",
    # variant nucleotide string, since all the dictionary keys
    # above are only (prefix, suffix) pairs which are missing these
    # nucleotides
    "variant_nucleotides",
]

VariantSequences = namedtuple(
    "VariantSequences",
    variant_sequences_object_fields)

def variant_reads_to_sequences(
        variant_reads,
        context_size=None,
        min_sequence_length=None,
        min_reads_per_sequence=DEFAULT_MIN_READS):
    """
    Parameters
    ----------
    variant_reads : list of VariantRead objects
        Each variant read should have the following fields:
            - prefix : str
            - variant : str
            - suffix : str
            - name : str

    context_size : int, optional
        Number of nucleotides left and right of the variant. If omitted then
        use full length of reads.

    min_sequence_length : int, optional
        Minimum length of detected sequence

    min_reads_per_sequence : int, optional
        Drop sequences which don't at least have this number of fully spanning
        reads.

    Returns a VariantSequences object, containing several different
    dictionaries mapping (prefix, variant, suffix) sequence tuples
    to names, counts, and weights of read support.
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
    logger.info("Sequence groups: %s" % (unique_sequence_groups,))

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
        if (len(prefix) + len(suffix) + variant_len) >= min_sequence_length and
        len(read_names) >= min_reads_per_sequence
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

    partially_supporting_read_names = defaultdict(set)
    partially_supporting_read_counts = Counter()
    partially_supporting_read_weights = Counter()

    for (prefix, suffix), reads in full_sequences.items():
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
                        assert fraction_bases_covered < 1.0, \
                            "Fraction for partial match must be <1.0, got %f " % (
                                fraction_bases_covered,)
                        partially_supporting_read_counts[other_key] += n_reads
                        partially_supporting_read_weights[other_key] += (
                            n_reads * fraction_bases_covered)
                        partially_supporting_read_names[other_key].update(reads)

    combined_weights = {
        key: full_count + partially_supporting_read_weights[key]
        for (key, full_count)
        in full_sequence_counts.items()
    }
    return VariantSequences(
        combined_sequence_weights=combined_weights,
        full_read_counts=full_sequence_counts,
        full_read_names=full_sequences,
        partial_read_counts=partially_supporting_read_counts,
        partial_read_weights=partially_supporting_read_weights,
        partial_read_names=partially_supporting_read_names,
        variant_nucleotides=variant_seq)

def variant_sequences_generator(
        variants,
        samfile,
        sequence_length=DEFAULT_SEQUENCE_LENGTH,
        min_reads=DEFAULT_MIN_READS):
    """
    For each variant, collect all possible sequence contexts around the
    variant which are spanned by at least min_reads.

    Parameters
    ----------
    variants : varcode.VariantCollection
        Variants for which we're trying to construct context sequences

    samfile : pysam.AlignmentFile
        BAM or SAM file containing RNA reads

    sequence_length : int
        Desired sequence length, including variant nucleotides

    min_reads : int
        Minimum number of reads supporting a particular sequence

    Generator that yields pairs of variants and VariantSequences objects.
    """
    for variant, variant_reads in variant_reads_generator(
            variants=variants,
            samfile=samfile):

        if len(variant_reads) == 0:
            logger.info("No variant reads found for %s" % variant)
            continue

        # the number of context nucleotides on either side of the variant
        # is half the desired length (minus the number of variant nucleotides)
        n_surrounding_nucleotides = sequence_length - len(variant.alt)

        flanking_context_size = int(np.ceil(n_surrounding_nucleotides / 2.0))
        logger.info(
            "Looking at %dnt RNA sequence context around %s" % (
                flanking_context_size,
                variant))

        sequences = variant_reads_to_sequences(
            variant_reads,
            context_size=flanking_context_size,
            min_reads_per_sequence=min_reads)

        yield variant, sequences

def variant_sequences_dataframe(
        variants,
        samfile,
        sequence_length=DEFAULT_SEQUENCE_LENGTH,
        min_reads=DEFAULT_MIN_READS):
    """
    Creates a dataframe of all detected cDNA sequences for the given variant
    collection and alignment file.

    Parameters
    ----------
    variants : varcode.VariantCollection
        Look for sequences containing these variants

    samfile : pysam.AlignmentFile
        Reads from which surrounding sequences are detected

    sequence_length : int
        Desired sequence length, including variant nucleotides

    min_reads : int
        Minimum number of reads supporting

    Returns pandas.DataFrame
    """
    columns = OrderedDict([
        ("chr", []),
        ("pos", []),
        ("ref", []),
        ("alt", [])
    ])
    for field in variant_sequences_object_fields:
        columns[field] = []

    for variant, sequences in variant_sequences_generator(
            variants=variants,
            samfile=samfile,
            sequence_length=sequence_length,
            min_reads=min_reads):
        for sequences_obj in sequences:
            columns["chr"].append(variant.contig)
            columns["pos"].append(variant.original_start)
            columns["ref"].append(variant.original_ref)
            columns["alt"].append(variant.original_alt)
            for field_name in variant_sequences_object_fields:
                columns[field_name] = getattr(sequences_obj, field_name)

    return pd.DataFrame(columns)
