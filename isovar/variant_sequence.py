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
from collections import namedtuple, OrderedDict
import logging

import numpy as np
import pandas as pd

from .common import group_unique_sequences, get_variant_nucleotides
from .variant_read import variant_reads_generator
from .default_parameters import (
    MIN_READS_SUPPORTING_VARIANT_CDNA_SEQUENCE,
    VARIANT_CDNA_SEQUENCE_LENGTH
)

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
        context_size=None,
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

    context_size : int, optional
        Number of nucleotides left and right of the variant. If omitted then
        use full length of reads.

    min_sequence_length : int, optional
        Minimum length of detected sequence

    min_reads_per_sequence : int, optional
        Drop sequences which don't at least have this number of fully spanning
        reads.

    Returns a collection of VariantSequence objects
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
        min_reads=MIN_READS_SUPPORTING_VARIANT_CDNA_SEQUENCE):
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

    Generator that yields pairs of variants and sorted list of VariantSequence
    objects.
    """
    for variant, variant_reads in variant_reads_generator(
            variants=variants,
            samfile=samfile):

        if len(variant_reads) == 0:
            logging.info("No variant reads found for %s" % variant)
            continue

        # the number of context nucleotides on either side of the variant
        # is half the desired length (minus the number of variant nucleotides)
        n_surrounding_nucleotides = sequence_length - len(variant.alt)

        flanking_context_size = int(np.ceil(n_surrounding_nucleotides / 2.0))
        logging.info(
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
        sequence_length=VARIANT_CDNA_SEQUENCE_LENGTH,
        min_reads=MIN_READS_SUPPORTING_VARIANT_CDNA_SEQUENCE):
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
    for field in VariantSequence._fields:
        columns[field] = []

    for variant, variant_sequences in variant_sequences_generator(
            variants=variants,
            samfile=samfile,
            sequence_length=sequence_length,
            min_reads=min_reads):
        for variant_sequence in variant_sequences:
            columns["chr"].append(variant.contig)
            columns["pos"].append(variant.original_start)
            columns["ref"].append(variant.original_ref)
            columns["alt"].append(variant.original_alt)
            for field_name in VariantSequence._fields:
                columns[field_name] = getattr(variant_sequence, field_name)

    return pd.DataFrame(columns)
