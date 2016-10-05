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

from .read_helpers import (
    get_single_allele_from_reads,
    group_unique_sequences,
)
from .variant_reads import filter_non_alt_reads_for_variant
from .default_parameters import (
    MIN_READS_SUPPORTING_VARIANT_CDNA_SEQUENCE,
    VARIANT_CDNA_SEQUENCE_LENGTH,
)
from .dataframe_builder import dataframe_from_generator


logger = logging.getLogger(__name__)


# using this base class to define the core fields of a VariantSequence
# but inheriting it from it to allow the addition of helper methods
VariantSequenceFields = namedtuple(
    "VariantSequence",
    [
        # nucleotides before a variant
        "prefix",
        # nucleotide sequence of a variant
        "alt",
        # nucleotides after a variant
        "suffix",
        # since we often want to look at prefix+alt+suffix, let's cache it
        "sequence",
        # reads which were used to determine this sequences
        "reads",
    ])

class VariantSequence(VariantSequenceFields):
    def __new__(cls, prefix, alt, suffix, reads):
        # construct sequence from prefix + alt + suffix
        return VariantSequenceFields.__new__(
            cls,
            prefix=prefix,
            alt=alt,
            suffix=suffix,
            sequence=prefix + alt + suffix,
            reads=frozenset(reads))

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
        The alleles must match and they must share at least `min_overlap_size`
        nucleotides in their prefix or suffix.
        """

        if self.alt != other.alt:
            # allele must match!
            return False

        if len(other.prefix) > len(self.prefix):
            # only consider strings that overlap like:
            #   variant_sequence1: ppppAssss
            #   variant_sequence2:   ppAsssssss
            # which excludes cases where variant_sequence2 has a longer
            # prefix
            return False
        elif len(other.suffix) < len(self.suffix):
            # similarly, we throw cases where variant_sequence2 is shorter
            # after the alt nucleotides than variant_sequence1
            return False

        # is the candidate sequence is a prefix of the accepted?
        # Example:
        # p1 a1 s1 = XXXXXXXX Y ZZZZZZ
        # p2 a2 s2 =       XX Y ZZZZZZZZZ
        # ...
        # then we can combine them into a longer sequence
        return (
            self.prefix.endswith(other.prefix) and
            other.suffix.startswith(self.suffix)
        )

    def add_reads(self, reads):
        """
        Create another VariantSequence with more supporting reads.
        """
        return VariantSequence(
            prefix=self.prefix,
            alt=self.alt,
            suffix=self.suffix,
            reads=self.reads.union(reads))

def sort_key_decreasing_read_count(variant_sequence):
    return -len(variant_sequence.reads)

def all_variant_sequences_supported_by_variant_reads(
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
        min_reads_supporting_cdna_sequence):
    n_total = len(variant_sequences)
    variant_sequences = [
        s
        for s in variant_sequences
        if len(s.reads) >= min_reads_supporting_cdna_sequence
    ]
    n_dropped = n_total - len(variant_sequences)
    if n_dropped > 0:
        logger.info("Dropped %d/%d variant sequences less than %d supporting reads",
            n_dropped,
            n_total,
            min_reads_supporting_cdna_sequence)
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
        logger.info("Dropped %d/%d variant sequences shorter than %d",
            n_dropped,
            n_total,
            min_required_sequence_length)
    return variant_sequences

def filter_variant_sequences(
        variant_sequences,
        preferred_sequence_length,
        min_reads_supporting_cdna_sequence=MIN_READS_SUPPORTING_VARIANT_CDNA_SEQUENCE):
    """
    Drop variant sequences which are shorter than request or don't have
    enough supporting reads.
    """
    variant_sequences = filter_variant_sequences_by_read_support(
        variant_sequences=variant_sequences,
        min_reads_supporting_cdna_sequence=min_reads_supporting_cdna_sequence)
    return filter_variant_sequences_by_length(
        variant_sequences=variant_sequences,
        preferred_sequence_length=preferred_sequence_length)

def supporting_reads_to_variant_sequences(
        variant_reads,
        preferred_sequence_length,
        min_reads_supporting_cdna_sequence=MIN_READS_SUPPORTING_VARIANT_CDNA_SEQUENCE):
    """
    Collapse variant-support RNA reads into consensus sequences of approximately
    the preferred length (may differ at the ends of transcripts), filter
    consensus sequences by length and number of supporting RNA reads.

    Parameters
    ----------
    variant_reads : list of AlleleRead objects
        Should all support the same variant allele nucleotides.

    preferred_sequence_length : int
        Total number of nucleotides in the assembled sequences, including
        variant nucleotides.

    min_sequence_length : int, optional
        Drop sequences shorter than this.

    min_reads_per_sequence : int, optional
        Drop sequences which don't at least have this number of fully spanning
        reads.

    Returns a collection of VariantSequence objects
    """
    # just in case variant_reads is a generator, convert it to a list
    variant_reads = list(variant_reads)

    if len(variant_reads) == 0:
        return []

    alt_seq = get_single_allele_from_reads(variant_reads)

    # the number of context nucleotides on either side of the variant
    # is half the desired length (minus the number of variant nucleotides)
    n_alt = len(alt_seq)
    n_surrounding_nucleotides = preferred_sequence_length - n_alt
    max_nucleotides_after_variant = n_surrounding_nucleotides // 2

    # if the number of nucleotides we need isn't divisible by 2 then
    # prefer to have one more *before* the variant since we need the
    # prefix sequence to match against reference transcripts
    max_nucleotides_before_variant = (
        n_surrounding_nucleotides - max_nucleotides_after_variant)

    variant_sequences = all_variant_sequences_supported_by_variant_reads(
        variant_reads=variant_reads,
        max_nucleotides_before_variant=max_nucleotides_before_variant,
        max_nucleotides_after_variant=max_nucleotides_after_variant)

    variant_sequences = filter_variant_sequences(
        variant_sequences=variant_sequences,
        preferred_sequence_length=preferred_sequence_length,
        min_reads_supporting_cdna_sequence=min_reads_supporting_cdna_sequence)

    # sort VariantSequence objects by decreasing order of supporting read
    # counts
    variant_sequences.sort(key=sort_key_decreasing_read_count)
    return variant_sequences

def overlapping_reads_to_variant_sequences(
        variant,
        overlapping_reads,
        min_reads_supporting_cdna_sequence,
        preferred_sequence_length):
    variant_reads = filter_non_alt_reads_for_variant(variant, overlapping_reads)
    return supporting_reads_to_variant_sequences(
        variant_reads,
        preferred_sequence_length=preferred_sequence_length,
        min_reads_supporting_cdna_sequence=min_reads_supporting_cdna_sequence)

def reads_generator_to_sequences_generator(
        variant_and_reads_generator,
        min_reads_supporting_cdna_sequence=MIN_READS_SUPPORTING_VARIANT_CDNA_SEQUENCE,
        preferred_sequence_length=VARIANT_CDNA_SEQUENCE_LENGTH):
    """
    For each variant, collect all possible sequence contexts around the
    variant which are spanned by at least min_reads.

    Parameters
    ----------
    variant_and_reads_generator : generator
        Sequence of Variant objects paired with a list of reads which
        overlap that variant.

    min_reads : int
        Minimum number of reads supporting variant sequence

    sequence_length : int
        Desired sequence length, including variant nucleotides

    Yields pairs with the following fields:
        - Variant
        - list of VariantSequence objects
    """
    for variant, variant_reads in variant_and_reads_generator:
        variant_sequences = overlapping_reads_to_variant_sequences(
            variant=variant,
            overlapping_reads=variant_reads,
            min_reads_supporting_cdna_sequence=min_reads_supporting_cdna_sequence,
            preferred_sequence_length=preferred_sequence_length)
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
