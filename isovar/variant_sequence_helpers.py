# Copyright (c) 2016-2019. Mount Sinai School of Medicine
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

"""
Helper functions for constructing and filtering VariantSequence objects
from reads overlapping a variant locus.
"""

from __future__ import print_function, division, absolute_import

from .allele_read_helpers import (
    get_single_allele_from_reads,
    group_unique_sequences,
    filter_non_alt_reads_for_variant
)
from .assembly import iterative_overlap_assembly, collapse_substrings
from .default_parameters import (
    MIN_VARIANT_SEQUENCE_COVERAGE,
    VARIANT_SEQUENCE_LENGTH,
    VARIANT_SEQUENCE_ASSEMBLY,
    MIN_VARIANT_SEQUENCE_ASSEMBLY_OVERLAP_SIZE,
)
from .variant_sequence import VariantSequence


def initial_variant_sequences_from_reads(
        variant_reads,
        max_nucleotides_before_variant=None,
        max_nucleotides_after_variant=None):
    """
    Get all unique sequences from reads spanning a variant locus. This will
    include partial sequences due to reads starting in the middle of the
    sequence around around a variant.
    """
    unique_sequence_groups = group_unique_sequences(
        variant_reads,
        max_prefix_size=max_nucleotides_before_variant,
        max_suffix_size=max_nucleotides_after_variant)

    return [
        VariantSequence(
            prefix=prefix,
            alt=alt,
            suffix=suffix,
            reads=reads)
        for ((prefix, alt, suffix), reads)
        in unique_sequence_groups.items()
    ]


def filter_variant_sequences_by_read_support(
        variant_sequences,
        min_variant_sequence_coverage):
    n_total = len(variant_sequences)
    variant_sequences = [
        s
        for s in variant_sequences
        if s.min_coverage() >= min_variant_sequence_coverage
    ]
    n_dropped = n_total - len(variant_sequences)
    if n_dropped > 0:
        logger.info(
            "Dropped %d/%d variant sequences less than %d supporting reads",
            n_dropped,
            n_total,
            min_variant_sequence_coverage)
    return variant_sequences


def filter_variant_sequences_by_length(
        variant_sequences,
        preferred_sequence_length):
    n_total = len(variant_sequences)
    if n_total == 0:
        return []
    # since we might have gotten some shorter fragments,
    # keep only the longest spanning sequence
    max_observed_sequence_length = max(len(s) for s in variant_sequences)

    # if we get back a sequence that's longer than the preferred length
    # then that doesn't mean we should necessarily drop the other sequences
    min_required_sequence_length = min(
        max_observed_sequence_length,
        preferred_sequence_length)

    variant_sequences = [
        s for s in variant_sequences
        if len(s.sequence) >= min_required_sequence_length
    ]
    n_dropped = n_total - len(variant_sequences)
    if n_dropped > 0:
        logger.info(
            "Dropped %d/%d variant sequences shorter than %d",
            n_dropped,
            n_total,
            min_required_sequence_length)
    return variant_sequences


def trim_variant_sequences(variant_sequences, min_variant_sequence_coverage):
    """
    Trim VariantSequences to desired coverage and then combine any
    subsequences which get generated.
    """
    n_total = len(variant_sequences)
    trimmed_variant_sequences = [
        variant_sequence.trim_by_coverage(min_variant_sequence_coverage)
        for variant_sequence in variant_sequences
    ]
    collapsed_variant_sequences = collapse_substrings(trimmed_variant_sequences)
    n_after_trimming = len(collapsed_variant_sequences)
    logger.info(
        "Kept %d/%d variant sequences after read coverage trimming to >=%dx",
        n_after_trimming,
        n_total,
        min_variant_sequence_coverage)
    return collapsed_variant_sequences


def filter_variant_sequences(
        variant_sequences,
        preferred_sequence_length,
        min_variant_sequence_coverage=MIN_VARIANT_SEQUENCE_COVERAGE,):
    """
    Drop variant sequences which are shorter than request or don't have
    enough supporting reads.
    """
    variant_sequences = trim_variant_sequences(
        variant_sequences, min_variant_sequence_coverage)

    return filter_variant_sequences_by_length(
        variant_sequences=variant_sequences,
        preferred_sequence_length=preferred_sequence_length)


def reads_to_variant_sequences(
        variant,
        reads,
        preferred_sequence_length,
        min_variant_sequence_coverage=MIN_VARIANT_SEQUENCE_COVERAGE,
        variant_sequence_assembly=VARIANT_SEQUENCE_ASSEMBLY):
    """
    Collapse variant-supporting RNA reads into consensus sequences of
    approximately the preferred length (may differ at the ends of transcripts),
    filter consensus sequences by length and number of supporting RNA reads.

    Parameters
    ----------
    variant : varcode.Variant

    reads : list of AlleleRead objects
        Should all support the same variant allele nucleotides.

    preferred_sequence_length : int
        Total number of nucleotides in the assembled sequences, including
        variant nucleotides.

    min_variant_sequence_coverage : int
        Drop sequences which don't at least have this number of reads
        covering each cDNA position.

    variant_sequence_assembly : bool
        Construct variant sequences by merging overlapping reads. If False
        then variant sequences must be fully spanned by cDNA reads.

    Returns a collection of VariantSequence objects
    """
    # just in case variant_reads is a generator, convert it to a list
    variant_reads = list(filter_non_alt_reads_for_variant(variant, reads))

    if len(variant_reads) == 0:
        return []

    alt_seq = get_single_allele_from_reads(variant_reads)

    # the number of context nucleotides on either side of the variant
    # is half the desired length (minus the number of variant nucleotides)
    n_alt_nucleotides = len(alt_seq)

    n_surrounding_nucleotides = preferred_sequence_length - n_alt_nucleotides
    max_nucleotides_after_variant = n_surrounding_nucleotides // 2

    # if the number of nucleotides we need isn't divisible by 2 then
    # prefer to have one more *before* the variant since we need the
    # prefix sequence to match against reference transcripts
    max_nucleotides_before_variant = (
        n_surrounding_nucleotides - max_nucleotides_after_variant)

    variant_sequences = initial_variant_sequences_from_reads(
        variant_reads=variant_reads,
        max_nucleotides_before_variant=max_nucleotides_before_variant,
        max_nucleotides_after_variant=max_nucleotides_after_variant)

    logger.info(
        "Initial pool of %d variant sequences (min length=%d, max length=%d)",
        len(variant_sequences),
        min(len(s) for s in variant_sequences),
        max(len(s) for s in variant_sequences))

    if variant_sequence_assembly:
        # this is a tricky parameter to set correctly:
        # by how many bases should two sequences overlap before
        # we merge, currently defaulting to either half the non-variant
        # nucleotides or 30 (whichever is smaller)
        variant_sequences = iterative_overlap_assembly(
            variant_sequences,
            min_overlap_size=min(
                MIN_VARIANT_SEQUENCE_ASSEMBLY_OVERLAP_SIZE,
                n_surrounding_nucleotides // 2))

    if variant_sequences:
        logger.info(
            "After overlap assembly: %d variant sequences (min length=%d, max length=%d)",
            len(variant_sequences),
            min(len(s) for s in variant_sequences),
            max(len(s) for s in variant_sequences))
    else:
        logger.info("After overlap assembly: 0 variant sequences")
        return []

    variant_sequences = filter_variant_sequences(
        variant_sequences=variant_sequences,
        preferred_sequence_length=preferred_sequence_length,
        min_variant_sequence_coverage=min_variant_sequence_coverage)

    if variant_sequences:
        logger.info(
            ("After coverage & length filtering: %d variant sequences "
             "(min length=%d, max length=%d)"),
            len(variant_sequences),
            min(len(s) for s in variant_sequences),
            max(len(s) for s in variant_sequences))
    else:
        logger.info("After coverage & length filtering: 0 variant sequences")
        return []

    # sort VariantSequence objects by decreasing order of supporting read
    # counts
    variant_sequences.sort(key=lambda vs: -len(vs.reads))
    return variant_sequences


def reads_generator_to_sequences_generator(
        variant_and_reads_generator,
        min_variant_sequence_coverage=MIN_VARIANT_SEQUENCE_COVERAGE,
        preferred_sequence_length=VARIANT_SEQUENCE_LENGTH,
        variant_sequence_assembly=VARIANT_SEQUENCE_ASSEMBLY):
    """
    For each variant, collect all possible sequence contexts around the
    variant which are spanned by at least min_reads.

    Parameters
    ----------
    variant_and_reads_generator : generator
        Sequence of Variant objects paired with a list of reads which
        overlap that variant.

    min_variant_sequence_coverage : int
        Minimum number of RNA reads supporting each nucleotide of the
        variant cDNA sequence

    preferred_sequence_length : int
        Desired sequence length, including variant nucleotides

    variant_sequence_assembly : bool
        Construct variant sequences by merging overlapping reads. If False
        then variant sequences must be fully spanned by cDNA reads.

    Yields pairs with the following fields:
        - Variant
        - list of VariantSequence objects
    """
    for variant, variant_reads in variant_and_reads_generator:
        variant_sequences = reads_to_variant_sequences(
            variant=variant,
            reads=variant_reads,
            min_variant_sequence_coverage=min_variant_sequence_coverage,
            preferred_sequence_length=preferred_sequence_length,
            variant_sequence_assembly=variant_sequence_assembly)
        yield variant, variant_sequences
