# Copyright (c) 2016-2018. Mount Sinai School of Medicine
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
This module combines variant cDNA sequences collected from a BAM file with
the reading frames of annotated reference transcripts to create candidate
translations.
"""

from __future__ import print_function, division, absolute_import

from six.moves import range, zip

from .dna import reverse_complement_dna
from .value_object import ValueObject
from .logging import get_logger

logger = get_logger(__name__)


class VariantSequenceInReadingFrame(ValueObject):
    """
    A variant cDNA sequence (possibly trimmed to get rid of low coverage
    tails) assigned to a particular strand ('+' or '-') and reading frame
    (one of 0, +1, +2). The strand and reading frame are determined by
    matching the cDNA sequence to a ReferenceContext.
    """
    __slots__ = [
        # since the reference context and variant sequence may have
        # different numbers of nucleotides before the variant, the cDNA prefix
        # gets truncated to the shortest length. To avoid having to recompute
        # that sequence again, let's just cache the full cDNA sequence we used
        # for translation here, along with an interval indicating which
        # nucleotides are from the variant of interest
        "cdna_sequence",
        "offset_to_first_complete_codon",
        "variant_cdna_interval_start",
        "variant_cdna_interval_end",
        "reference_cdna_sequence_before_variant",
        "reference_cdna_sequence_after_variant",
        "number_mismatches_before_variant",
        "number_mismatches_after_variant"
    ]

    def __init__(
            self,
            cdna_sequence,
            offset_to_first_complete_codon,
            variant_cdna_interval_start,
            variant_cdna_interval_end,
            reference_cdna_sequence_before_variant,
            reference_cdna_sequence_after_variant,
            number_mismatches_before_variant,
            number_mismatches_after_variant):
        self.cdna_sequence = cdna_sequence
        self.offset_to_first_complete_codon = offset_to_first_complete_codon
        self.variant_cdna_interval_start = variant_cdna_interval_start
        self.variant_cdna_interval_end = variant_cdna_interval_end
        self.reference_cdna_sequence_before_variant = (
            reference_cdna_sequence_before_variant)
        self.reference_cdna_sequence_after_variant = (
            reference_cdna_sequence_after_variant)
        self.number_mismatches_before_variant = number_mismatches_before_variant
        self.number_mismatches_after_variant = number_mismatches_after_variant

    @property
    def in_frame_cdna_sequence(self):
        from_first_codon = self.cdna_sequence[self.offset_to_first_complete_codon:]
        # get rid of any trailing out of frame nucleotides
        in_frame_length = len(from_first_codon) // 3 * 3
        return from_first_codon[:in_frame_length]

    @classmethod
    def from_variant_sequence_and_reference_context(
            cls, variant_sequence, reference_context):
        """
        Combines a VariantSequence with the reading frame implied by a
        ReferenceContext, reverse complementing if necessary and finding the
        offset to the first complete codon in the cDNA sequence.
        Parameters
        ----------
        variant_sequence : VariantSequence

        reference_context : ReferenceContext

        Returns a VariantSequenceInReadingFrame object
        """
        (cdna_prefix, cdna_alt, cdna_suffix,
            reference_prefix, reference_suffix, n_trimmed_from_reference) = trim_sequences(
                variant_sequence, reference_context)

        logger.info(
            ("cdna_predix='%s', cdna_alt='%s', cdna_suffix='%s', "
             "reference_prefix='%s', reference_suffix='%s', n_trimmed=%d"),
            cdna_prefix,
            cdna_alt,
            cdna_suffix,
            reference_prefix,
            reference_suffix,
            n_trimmed_from_reference)

        n_mismatch_before_variant = count_mismatches_before_variant(reference_prefix, cdna_prefix)
        n_mismatch_after_variant = count_mismatches_after_variant(reference_suffix, cdna_suffix)

        ref_codon_offset = reference_context.offset_to_first_complete_codon

        # ReferenceContext carries with an offset to the first complete codon
        # in the reference sequence. This may need to be adjusted if the reference
        # sequence is longer than the variant sequence (and thus needs to be trimmed)
        offset_to_first_complete_codon = compute_offset_to_first_complete_codon(
            offset_to_first_complete_reference_codon=ref_codon_offset,
            n_trimmed_from_reference_sequence=n_trimmed_from_reference)

        cdna_sequence = cdna_prefix + cdna_alt + cdna_suffix
        variant_interval_start = len(cdna_prefix)
        variant_interval_end = variant_interval_start + len(cdna_alt)

        return VariantSequenceInReadingFrame(
            cdna_sequence=cdna_sequence,
            offset_to_first_complete_codon=offset_to_first_complete_codon,
            variant_cdna_interval_start=variant_interval_start,
            variant_cdna_interval_end=variant_interval_end,
            reference_cdna_sequence_before_variant=reference_prefix,
            reference_cdna_sequence_after_variant=reference_suffix,
            number_mismatches_before_variant=n_mismatch_before_variant,
            number_mismatches_after_variant=n_mismatch_after_variant)


def trim_sequences(variant_sequence, reference_context):
    """
    A VariantSequence and ReferenceContext may contain a different number of
    nucleotides before the variant locus. Furthermore, the VariantSequence is
    always expressed in terms of the positive strand against which it aligned,
    but reference transcripts may have sequences from the negative strand of the
    genome. Take the reverse complement of the VariantSequence if the
    ReferenceContext is from negative strand transcripts and trim either
    sequence to ensure that the prefixes are of the same length.

    Parameters
    ----------
    variant_sequence : VariantSequence

    reference_context : ReferenceContext

    Returns a tuple with the following fields:
        1) cDNA prefix of variant sequence, trimmed to be same length as the
           reference prefix. If the reference context was on the negative
           strand then this is the trimmed sequence *after* the variant from
           the genomic DNA sequence.
        2) cDNA sequence of the variant nucleotides, in reverse complement if
           the reference context is on the negative strand.
        3) cDNA sequence of the nucleotides after the variant nucleotides. If
           the reference context is on the negative strand then this sequence
           is the reverse complement of the original prefix sequence.
        4) Reference sequence before the variant locus, trimmed to be the
           same length as the variant prefix.
        5) Reference sequence after the variant locus, untrimmed.
        6) Number of nucleotides trimmed from the reference sequence, used
           later for adjustint offset to first complete codon.
    """
    cdna_prefix = variant_sequence.prefix
    cdna_alt = variant_sequence.alt
    cdna_suffix = variant_sequence.suffix

    # if the transcript is on the reverse strand then we have to
    # take the sequence PREFIX|VARIANT|SUFFIX
    # and take the complement of XIFFUS|TNAIRAV|XIFERP
    if reference_context.strand == "-":
        # notice that we are setting the *prefix* to be reverse complement
        # of the *suffix* and vice versa
        cdna_prefix, cdna_alt, cdna_suffix = (
            reverse_complement_dna(cdna_suffix),
            reverse_complement_dna(cdna_alt),
            reverse_complement_dna(cdna_prefix)
        )

    reference_sequence_before_variant = reference_context.sequence_before_variant_locus
    reference_sequence_after_variant = reference_context.sequence_after_variant_locus

    # trim the reference prefix and the RNA-derived prefix sequences to the same length
    if len(reference_sequence_before_variant) > len(cdna_prefix):
        n_trimmed_from_reference = len(reference_sequence_before_variant) - len(cdna_prefix)
        n_trimmed_from_variant = 0
    elif len(reference_sequence_before_variant) < len(cdna_prefix):
        n_trimmed_from_variant = len(cdna_prefix) - len(reference_sequence_before_variant)
        n_trimmed_from_reference = 0
    else:
        n_trimmed_from_variant = 0
        n_trimmed_from_reference = 0

    reference_sequence_before_variant = reference_sequence_before_variant[
        n_trimmed_from_reference:]
    cdna_prefix = cdna_prefix[n_trimmed_from_variant:]

    return (
        cdna_prefix,
        cdna_alt,
        cdna_suffix,
        reference_sequence_before_variant,
        reference_sequence_after_variant,
        n_trimmed_from_reference
    )


def count_mismatches_before_variant(reference_prefix, cdna_prefix):
    """
    Computes the number of mismatching nucleotides between two cDNA sequences before a variant
    locus.

    Parameters
    ----------
    reference_prefix : str
        cDNA sequence of a reference transcript before a variant locus

    cdna_prefix : str
        cDNA sequence detected from RNAseq before a variant locus
    """
    if len(reference_prefix) != len(cdna_prefix):
        raise ValueError(
            "Expected reference prefix '%s' to be same length as %s" % (
                reference_prefix, cdna_prefix))
    return sum(xi != yi for (xi, yi) in zip(reference_prefix, cdna_prefix))


def count_mismatches_after_variant(reference_suffix, cdna_suffix):
    """
    Computes the number of mismatching nucleotides between two cDNA sequences after a variant locus.

    Parameters
    ----------
    reference_suffix : str
        cDNA sequence of a reference transcript after a variant locus

    cdna_suffix : str
        cDNA sequence detected from RNAseq after a variant locus
    """

    len_diff = len(cdna_suffix) - len(reference_suffix)

    # if the reference is shorter than the read, the read runs into the intron - these count as
    # mismatches
    return sum(xi != yi for (xi, yi) in zip(reference_suffix, cdna_suffix)) + max(0, len_diff)


def compute_offset_to_first_complete_codon(
        offset_to_first_complete_reference_codon,
        n_trimmed_from_reference_sequence):
    """
    Once we've aligned the variant sequence to the ReferenceContext, we need
    to transfer reading frame from the reference transcripts to the variant
    sequences.

    Parameters
    ----------
    offset_to_first_complete_reference_codon : int

    n_trimmed_from_reference_sequence : int

    Returns an offset into the variant sequence that starts from a complete
    codon.
    """
    if n_trimmed_from_reference_sequence <= offset_to_first_complete_reference_codon:
        return (
            offset_to_first_complete_reference_codon -
            n_trimmed_from_reference_sequence)
    else:
        n_nucleotides_trimmed_after_first_codon = (
            n_trimmed_from_reference_sequence -
            offset_to_first_complete_reference_codon)
        frame = n_nucleotides_trimmed_after_first_codon % 3
        return (3 - frame) % 3


def match_variant_sequence_to_reference_context(
        variant_sequence,
        reference_context,
        min_transcript_prefix_length,
        max_transcript_mismatches,
        include_mismatches_after_variant=False,
        max_trimming_attempts=2):
    """
    Iteratively trim low-coverage subsequences of a variant sequence
    until it either matches the given reference context or there
    are too few nucleotides left in the variant sequence.

    Parameters
    ----------
    variant_sequence : VariantSequence
        Assembled sequence from RNA reads, will need to be to be reverse
        complemented if matching against a reference transcript on the
        negative strand.

    reference_context : ReferenceContext
        Sequence of reference transcript before the variant and associated
        metadata.

    min_transcript_prefix_length : int
        Minimum number of nucleotides we try to match against a reference
        transcript.

    max_transcript_mismatches : int
        Maximum number of nucleotide differences between reference transcript
        sequence and the variant sequence.

    include_mismatches_after_variant : bool
        Set to true if the number of mismatches after the variant locus should
        count toward the total max_transcript_mismatches, which by default
        only counts mismatches before the variant locus.

    max_trimming_attempts : int
        How many times do we try trimming the VariantSequence to higher
        levels of coverage before giving up?

    Returns VariantSequenceInReadingFrame or None
    """
    variant_sequence_in_reading_frame = None

    # if we can't get the variant sequence to match this reference
    # context then keep trimming it by coverage until either
    for i in range(max_trimming_attempts + 1):
        # check the reverse-complemented prefix if the reference context is
        # on the negative strand since variant sequence is aligned to
        # genomic DNA (positive strand)
        variant_sequence_too_short = (
            (reference_context.strand == "+" and
                len(variant_sequence.prefix) < min_transcript_prefix_length) or
            (reference_context.strand == "-" and
                len(variant_sequence.suffix) < min_transcript_prefix_length)
        )
        if variant_sequence_too_short:
            logger.info(
                "Variant sequence %s shorter than min allowed %d (iter=%d)",
                variant_sequence,
                min_transcript_prefix_length,
                i + 1)
            return None

        variant_sequence_in_reading_frame = \
            VariantSequenceInReadingFrame.from_variant_sequence_and_reference_context(
                variant_sequence=variant_sequence,
                reference_context=reference_context)

        if variant_sequence_in_reading_frame is None:
            return None

        n_mismatch_before_variant = (
            variant_sequence_in_reading_frame.number_mismatches_before_variant)
        n_mismatch_after_variant = (
            variant_sequence_in_reading_frame.number_mismatches_after_variant)

        logger.info("Iter #%d/%d: %s" % (
            i + 1,
            max_trimming_attempts + 1,
            variant_sequence_in_reading_frame))

        total_mismatches = n_mismatch_before_variant
        if include_mismatches_after_variant:
            total_mismatches += n_mismatch_after_variant
        if total_mismatches <= max_transcript_mismatches:
            # if we got a variant sequence + reading frame with sufficiently
            # few mismatches then call it a day
            return variant_sequence_in_reading_frame

        logger.info(
            ("Too many mismatches (%d) between variant sequence %s and "
             "reference context %s (attempt=%d/%d)"),
            n_mismatch_before_variant,
            variant_sequence,
            reference_context,
            i + 1,
            max_trimming_attempts + 1)
        # if portions of the sequence are supported by only 1 read
        # then try trimming to 2 to see if the better supported
        # subsequence can be better matched against the reference
        current_min_coverage = variant_sequence.min_coverage()
        logger.info(
            "Trimming to subsequence covered by at least %d reads",
            current_min_coverage + 1)
        variant_sequence = variant_sequence.trim_by_coverage(
            current_min_coverage + 1)
    return None
