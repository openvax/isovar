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
from collections import namedtuple
import logging


from .common import group_unique_sequences, get_variant_nucleotides
from .variant_read import variant_reads_generator
from .default_parameters import (
    MIN_READS_SUPPORTING_VARIANT_CDNA_SEQUENCE,
    VARIANT_CDNA_SEQUENCE_LENGTH,
    MIN_READ_MAPPING_QUALITY,
)
from .dataframe_builder import DataFrameBuilder

VariantSequence = namedtuple(
    "VariantSequence",
    [
        "prefix",
        "alt",
        "suffix",
        # combined sequence cached for convenience, so we don't have to
        # repeatedly concatenate prefix + variant_nucleotides + suffix
        "full_sequence",
        "read_names",
        "read_count",
    ]
)

def sort_key_decreasing_read_count(variant_sequence):
    return -variant_sequence.read_count

def variant_reads_to_sequences(
        variant_reads,
        max_nucleotides_before_variant=None,
        max_nucleotides_after_variant=None,
        min_sequence_length=None,
        min_reads_per_sequence=MIN_READS_SUPPORTING_VARIANT_CDNA_SEQUENCE):
    """
    Parameters
    ----------
    variant_reads : list of VariantRead objects
        Each variant read should have the following fields:
            - prefix : str
            - variant : str
            - suffix : str
            - name : str

    max_nucleotides_before_variant : int, optional
        Number of nucleotides left of a variant. If omitted then
        use all prefix nucleotides.

    max_nucleotides_after_variant : int, optional
        Number of nucleotides right of a variant. If omitted then
        use all suffix nucleotides.

    min_sequence_length : int, optional
        Minimum length of detected sequence

    min_reads_per_sequence : int, optional
        Drop sequences which don't at least have this number of fully spanning
        reads.

    Returns a collection of VariantSequence objects
    """
    # just in case variant_reads is a generator, convert it to a list
    variant_reads = list(variant_reads)

    if len(variant_reads) == 0:
        return []

    # Get all unique sequences from reads spanning the
    # variant locus. This will include partial sequences
    # due to reads starting in the middle of the context,
    # we'll use these only to compute partial support for a
    # full length sequence
    unique_sequence_groups = group_unique_sequences(
        variant_reads,
        max_prefix_size=max_nucleotides_before_variant,
        max_suffix_size=max_nucleotides_after_variant)

    variant_seq = get_variant_nucleotides(variant_reads)
    variant_len = len(variant_seq)

    if min_sequence_length is None or min_sequence_length == 0:
        # if not specified, then use the longest sequence available
        # which will be at most twice the context size
        min_sequence_length = variant_len + max(
            len(prefix) + len(suffix)
            for (prefix, suffix) in unique_sequence_groups.keys())

    variant_sequences = [
        VariantSequence(
            prefix=prefix,
            alt=variant_seq,
            suffix=suffix,
            full_sequence=prefix + variant_seq + suffix,
            read_names=read_names,
            read_count=len(read_names))
        for ((prefix, suffix), read_names)
        in unique_sequence_groups.items()
    ]
    n_total = len(variant_sequences)

    variant_sequences = [
        x for x in variant_sequences
        if len(x.full_sequence) >= min_sequence_length and
        len(x.read_names) >= min_reads_per_sequence
    ]

    n_dropped = n_total - len(variant_sequences)
    logging.info("Dropped %d/%d variant sequences" % (n_dropped, n_total))

    # sort VariantSequence objects by decreasing order of supporting read
    # counts
    variant_sequences.sort(key=sort_key_decreasing_read_count)
    return variant_sequences

def variant_sequences_generator(
        variants,
        samfile,
        sequence_length=VARIANT_CDNA_SEQUENCE_LENGTH,
        min_reads=MIN_READS_SUPPORTING_VARIANT_CDNA_SEQUENCE,
        min_mapping_quality=MIN_READ_MAPPING_QUALITY):
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

    min_mapping_quality : int
        Minimum MAPQ value before a read gets ignored

    Generator that yields tuples with the following fields:
        - Variant
        - list of VariantSequence objects
    """
    for variant, variant_reads in variant_reads_generator(
            variants=variants,
            samfile=samfile,
            min_mapping_quality=min_mapping_quality):
        # the number of context nucleotides on either side of the variant
        # is half the desired length (minus the number of variant nucleotides)
        n_alt = len(variant.alt)
        n_surrounding_nucleotides = sequence_length - n_alt
        max_nucleotides_after_variant = n_surrounding_nucleotides // 2
        # if the number of nucleotides we need isn't divisible by 2 then
        # prefer to have one more *before* the variant since we need the
        # prefix sequence to match against reference transcripts
        max_nucleotides_before_variant = (
            n_surrounding_nucleotides - max_nucleotides_after_variant)
        logging.info(
            "Looking at %dnt before and %dnt after variant %s" % (
                max_nucleotides_before_variant,
                max_nucleotides_after_variant,
                variant))

        variant_sequences = variant_reads_to_sequences(
            variant_reads,
            max_nucleotides_before_variant=max_nucleotides_before_variant,
            max_nucleotides_after_variant=max_nucleotides_after_variant,
            min_reads_per_sequence=min_reads)

        yield variant, variant_sequences

def variant_sequences_dataframe(
        variants,
        samfile,
        sequence_length=VARIANT_CDNA_SEQUENCE_LENGTH,
        min_reads=MIN_READS_SUPPORTING_VARIANT_CDNA_SEQUENCE,
        min_mapping_quality=MIN_READ_MAPPING_QUALITY):
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

    min_mapping_quality : int
        Minimum MAPQ value before a read gets ignored

    Returns pandas.DataFrame
    """
    df_builder = DataFrameBuilder(VariantSequence)
    for variant, variant_sequences, variant_reads in variant_sequences_generator(
            variants=variants,
            samfile=samfile,
            sequence_length=sequence_length,
            min_reads=min_reads,
            min_mapping_quality=min_mapping_quality):
        for variant_sequence in variant_sequences:
            df_builder.add(variant, variant_sequence)
    return df_builder.to_dataframe()
