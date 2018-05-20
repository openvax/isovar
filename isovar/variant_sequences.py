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

import numpy as np

from .read_helpers import (
    get_single_allele_from_reads,
    group_unique_sequences,
)
from .variant_reads import filter_non_alt_reads_for_variant
from .default_parameters import (
    MIN_ALT_RNA_READS,
    MIN_VARIANT_SEQUENCE_COVERAGE,
    VARIANT_SEQUENCE_LENGTH,
    VARIANT_SEQUENCE_ASSEMBLY,
    MIN_VARIANT_SEQUENCE_ASSEMBLY_OVERLAP_SIZE,
)
from .dataframe_builder import dataframe_from_generator
from .assembly import iterative_overlap_assembly, collapse_substrings
from .value_object import ValueObject
from .logging import get_logger

logger = get_logger(__name__)


class VariantSequence(ValueObject):
    __slots__ = [
        # nucleotides before a variant
        "prefix",
        # nucleotide sequence of a variant
        "alt",
        # nucleotides after a variant
        "suffix",
        # since we often want to look at prefix+alt+suffix, let's cache it
        "sequence",
        # reads which were used to determine this sequences
        "reads"
    ]

    def __init__(self, prefix, alt, suffix, reads):
        self.prefix = prefix
        self.alt = alt
        self.suffix = suffix
        self.sequence = prefix + alt + suffix
        self.reads = frozenset(reads)

    def __len__(self):
        return len(self.sequence)

    @property
    def read_names(self):
        return {r.name for r in self.reads}

    def contains(self, other):
        """
        Is the other VariantSequence a subsequence of this one?

        The two sequences must agree on the alt nucleotides, the prefix of the
        longer must contain the prefix of the shorter, and the suffix of the
        longer must contain the suffix of the shorter.
        """
        return (self.alt == other.alt and
                self.prefix.endswith(other.prefix) and
                self.suffix.startswith(other.suffix))

    def left_overlaps(self, other, min_overlap_size=1):
        """
        Does this VariantSequence overlap another on the left side?
        """

        if self.alt != other.alt:
            # allele must match!
            return False

        if len(other.prefix) > len(self.prefix):
            # only consider strings that overlap like:
            #   self: ppppAssss
            #   other:   ppAsssssss
            # which excludes cases where the other sequence has a longer
            # prefix
            return False
        elif len(other.suffix) < len(self.suffix):
            # similarly, we throw away cases where the other sequence is shorter
            # after the alt nucleotides than this sequence
            return False

        # is the other sequence a prefix of this sequence?
        # Example:
        # p1 a1 s1 = XXXXXXXX Y ZZZZZZ
        # p2 a2 s2 =       XX Y ZZZZZZZZZ
        # ...
        # then we can combine them into a longer sequence
        sequence_overlaps = (
            self.prefix.endswith(other.prefix) and
            other.suffix.startswith(self.suffix)
        )
        prefix_overlap_size = min(len(self.prefix), len(other.prefix))
        suffix_overlap_size = min(len(other.suffix), len(self.suffix))
        overlap_size = (
            prefix_overlap_size + suffix_overlap_size + len(self.alt))

        return sequence_overlaps and overlap_size >= min_overlap_size

    def add_reads(self, reads):
        """
        Create another VariantSequence with more supporting reads.
        """
        if len(reads) == 0:
            return self
        new_reads = self.reads.union(reads)
        if len(new_reads) > len(self.reads):
            return VariantSequence(
                prefix=self.prefix,
                alt=self.alt,
                suffix=self.suffix,
                reads=new_reads)
        else:
            return self

    def combine(self, other_sequence):
        """
        If this sequence is the prefix of another sequence, combine
        them into a single VariantSequence object. If the other sequence
        is contained in this one, then add its reads to this VariantSequence.
        Also tries to flip the order (e.g. this sequence is a suffix or
        this sequence is a subsequence). If sequences can't be combined
        then returns None.
        """
        if other_sequence.alt != self.alt:
            logger.warn(
                "Cannot combine %s and %s with mismatching alt sequences",
                self,
                other_sequence)
            return None
        elif self.contains(other_sequence):
            return self.add_reads(other_sequence.reads)
        elif other_sequence.contains(self):
            return other_sequence.add_reads(self.reads)
        elif self.left_overlaps(other_sequence):
            # If sequences are like AABC and ABCC
            return VariantSequence(
                prefix=self.prefix,
                alt=self.alt,
                suffix=other_sequence.suffix,
                reads=self.reads.union(other_sequence.reads))
        elif other_sequence.left_overlaps(self):
            return VariantSequence(
                prefix=other_sequence.prefix,
                alt=self.alt,
                suffix=self.suffix,
                reads=self.reads.union(other_sequence.reads))
        else:
            # sequences don't overlap
            return None

    def variant_indices(self):
        """
        When we combine prefix + alt + suffix into a single string,
        what are is base-0 index interval which gets us back the alt
        sequence? First returned index is inclusive, the second is exclusive.
        """
        variant_start_index = len(self.prefix)
        variant_len = len(self.alt)
        variant_end_index = variant_start_index + variant_len
        return variant_start_index, variant_end_index

    def coverage(self):
        """
        Returns NumPy array indicating number of reads covering each
        nucleotides of this sequence.
        """
        variant_start_index, variant_end_index = self.variant_indices()
        n_nucleotides = len(self)
        coverage_array = np.zeros(n_nucleotides, dtype="int32")
        for read in self.reads:
            coverage_array[
                max(0, variant_start_index - len(read.prefix)):
                min(n_nucleotides, variant_end_index + len(read.suffix))] += 1
        return coverage_array

    def min_coverage(self):
        return np.min(self.coverage())

    def mean_coverage(self):
        return np.mean(self.coverage())

    def trim_by_coverage(self, min_reads):
        """
        Given the min number of reads overlapping each nucleotide of
        a variant sequence, trim this sequence by getting rid of positions
        which are overlapped by fewer reads than specified.
        """
        read_count_array = self.coverage()
        logger.info("Coverage: %s (len=%d)" % (
            read_count_array, len(read_count_array)))
        sufficient_coverage_mask = read_count_array >= min_reads
        sufficient_coverage_indices = np.argwhere(sufficient_coverage_mask)
        if len(sufficient_coverage_indices) == 0:
            logger.debug("No bases in %s have coverage >= %d" % (self, min_reads))
            return VariantSequence(prefix="", alt="", suffix="", reads=self.reads)
        variant_start_index, variant_end_index = self.variant_indices()
        # assuming that coverage drops off monotonically away from
        # variant nucleotides
        first_covered_index = sufficient_coverage_indices.min()
        last_covered_index = sufficient_coverage_indices.max()
        # adding 1 to last_covered_index since it's an inclusive index
        # whereas variant_end_index is the end of a half-open interval
        if (first_covered_index > variant_start_index or
                last_covered_index + 1 < variant_end_index):
            # Example:
            #   Nucleotide sequence:
            #       ACCCTTTT|AA|GGCGCGCC
            #   Coverage:
            #       12222333|44|33333211
            # Then the mask for bases covered >= 4x would be:
            #       ________|**|________
            # with indices:
            #       first_covered_index = 9
            #       last_covered_index = 10
            #       variant_start_index = 9
            #       variant_end_index = 11
            logger.debug("Some variant bases in %s don't have coverage >= %d" % (
                self, min_reads))
            return VariantSequence(prefix="", alt="", suffix="", reads=self.reads)
        return VariantSequence(
            prefix=self.prefix[first_covered_index:],
            alt=self.alt,
            suffix=self.suffix[:last_covered_index - variant_end_index + 1],
            reads=self.reads)


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
        min_alt_rna_reads=MIN_ALT_RNA_READS,
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

    min_alt_rna_reads : int
        Drop sequences from loci which lack at least this number of reads
        that agree with the variant allele.

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

    if len(variant_reads) < min_alt_rna_reads:
        logger.info(
            "Skipping %s because only %d alt RNA reads (min=%d)",
            variant,
            len(variant_reads),
            min_alt_rna_reads)
        return []
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
        min_alt_rna_reads=MIN_ALT_RNA_READS,
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

    min_alt_rna_reads : int
        Minimum number of RNA reads supporting variant allele

    min_variant_sequence_coverage : int
        Minimum number of RNA reads supporting each nucleotide of the
        variant cDNA sequence

    sequence_length : int
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
            min_alt_rna_reads=min_alt_rna_reads,
            min_variant_sequence_coverage=min_variant_sequence_coverage,
            preferred_sequence_length=preferred_sequence_length,
            variant_sequence_assembly=variant_sequence_assembly)
        yield variant, variant_sequences


def variant_sequences_generator_to_dataframe(variant_sequences_generator):
    """
    Creates a dataframe from a generator which yields
    (Variant, [VariantSequence]) pairs.

    Returns pandas.DataFrame
    """
    # TODO: Change VariantSequence.alt to VariantSequence.alt_nucleotides
    # or something else that doesn't clash with a variant's `alt` field
    return dataframe_from_generator(
        VariantSequence,
        variant_sequences_generator,
        rename_dict={"alt": "allele"},
        extra_column_fns={
            "gene": lambda variant, _: ";".join(variant.gene_names),
        })
