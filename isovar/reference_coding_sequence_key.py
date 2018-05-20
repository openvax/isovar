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

from .variant_helpers import interbase_range_affected_by_variant_on_transcript
from .reference_sequence_key import ReferenceSequenceKey
from .logging import get_logger

logger = get_logger(__name__)


class ReferenceCodingSequenceKey(ReferenceSequenceKey):
    """
    ReferenceCodingSequenceKey includes all the fields of a ReferenceSequenceKey,
    and additionally tracks the reading frame and information about where the
    start codon and 5' UTR are relative to this sequence fragment.
    """

    # additional fields on top of ReferenceSequenceKey
    __slots__ = [
        # if the reference context includes the 5' UTR then
        # this is the offset to the start codon, otherwise it's the
        # offset needed to get the first base of a codon
        "offset_to_first_complete_codon",
        # does this context overlap a start codon?
        "overlaps_start_codon",
        # does this context contain the whole trinucleotide start codon?
        "contains_start_codon",
        # does this context contain any UTR bases?
        "contains_five_prime_utr",
        # translation of complete codons in the reference context
        # before the variant
        "amino_acids_before_variant"
    ]

    def __init__(
            self,
            strand,
            sequence_before_variant_locus,
            sequence_at_variant_locus,
            sequence_after_variant_locus,
            offset_to_first_complete_codon,
            contains_start_codon,
            overlaps_start_codon,
            contains_five_prime_utr,
            amino_acids_before_variant):
        ReferenceSequenceKey.__init__(
            self,
            strand=strand,
            sequence_before_variant_locus=sequence_before_variant_locus,
            sequence_at_variant_locus=sequence_at_variant_locus,
            sequence_after_variant_locus=sequence_after_variant_locus)
        self.offset_to_first_complete_codon = offset_to_first_complete_codon
        self.overlaps_start_codon = overlaps_start_codon
        self.contains_start_codon = contains_start_codon
        self.contains_five_prime_utr = contains_five_prime_utr
        self.amino_acids_before_variant = amino_acids_before_variant

    @classmethod
    def from_variant_and_transcript_and_sequence_key(
            cls, variant, transcript, sequence_key):
        """
        Assuming that the transcript has a coding sequence, take a
        ReferenceSequenceKey (region of the transcript around the variant) and
        return a ReferenceCodingSequenceKey (or None).
        """
        # get the interbase range of offsets which capture all reference
        # bases modified by the variant
        variant_start_offset, variant_end_offset = \
            interbase_range_affected_by_variant_on_transcript(
                variant=variant,
                transcript=transcript)

        start_codon_idx = min(transcript.start_codon_spliced_offsets)

        # skip any variants which occur in the 5' UTR or overlap the start codon
        # since
        #   (1) UTR variants have unclear coding effects and
        #   (2) start-loss variants may result in a new start codon / reading
        #       frame but we're not sure which one!
        if variant_start_offset < start_codon_idx + 3:
            logger.info(
                "Skipping transcript %s for variant %s, must be after start codon",
                transcript.name,
                variant)
            return None

        stop_codon_idx = min(transcript.stop_codon_spliced_offsets)

        # skip variants which affect the 3' UTR of the transcript since
        # they don't have obvious coding effects on the protein sequence
        if variant_start_offset >= stop_codon_idx + 3:
            logger.info(
                "Skipping transcript %s for variant %s, occurs in 3' UTR",
                transcript,
                variant)
            return None

        n_prefix = len(sequence_key.sequence_before_variant_locus)
        prefix_start_idx = variant_start_offset - n_prefix
        n_bases_between_start_and_variant = variant_start_offset - start_codon_idx
        n_full_codons_before_variant = n_bases_between_start_and_variant // 3

        # if the sequence before the variant contains more bases than the
        # distance to the start codon, then by definition it must contain
        # some untranslated bases
        contains_five_prime_utr = (n_prefix > n_bases_between_start_and_variant)
        # allows for the possibility that the first base in the sequence might
        # be the first nucleotide of the start codon
        contains_start_codon = (n_prefix >= n_bases_between_start_and_variant)
        # the sequence context might only include the 2nd or 3rd bases of
        # the start codon
        overlaps_start_codon = (n_prefix > n_bases_between_start_and_variant - 3)

        if contains_start_codon:
            offset_to_first_complete_codon = start_codon_idx - prefix_start_idx
            amino_acids_before_variant = transcript.protein_sequence[:n_full_codons_before_variant]
        else:
            reading_frame = (prefix_start_idx - start_codon_idx) % 3
            offset_to_first_complete_codon = reading_frame_to_offset(reading_frame)
            n_codons_in_prefix = (n_prefix - offset_to_first_complete_codon) // 3
            amino_acids_before_variant = transcript.protein_sequence[
                n_full_codons_before_variant - n_codons_in_prefix:
                n_full_codons_before_variant]

        return cls(
            strand=sequence_key.strand,
            sequence_before_variant_locus=sequence_key.sequence_before_variant_locus,
            sequence_at_variant_locus=sequence_key.sequence_at_variant_locus,
            sequence_after_variant_locus=sequence_key.sequence_after_variant_locus,
            offset_to_first_complete_codon=offset_to_first_complete_codon,
            contains_start_codon=contains_start_codon,
            overlaps_start_codon=overlaps_start_codon,
            contains_five_prime_utr=contains_five_prime_utr,
            amino_acids_before_variant=amino_acids_before_variant)

    @classmethod
    def from_variant_and_transcript(
            cls,
            variant,
            transcript,
            context_size):
        """
        Extracts the reference sequence around a variant locus on a particular
        transcript and determines the reading frame at the start of that
        sequence context.

        Parameters
        ----------
        variant : varcode.Variant

        transcript : pyensembl.Transcript

        context_size : int

        Returns SequenceKeyWithReadingFrame object or None if Transcript lacks
        coding sequence, protein sequence or annotated start/stop codons.
        """
        if not transcript.contains_start_codon:
            logger.info(
                "Expected transcript %s for variant %s to have start codon",
                transcript.name,
                variant)
            return None

        if not transcript.contains_stop_codon:
            logger.info(
                "Expected transcript %s for variant %s to have stop codon",
                transcript.name,
                variant)
            return None

        if not transcript.protein_sequence:
            logger.info(
                "Expected transript %s for variant %s to have protein sequence",
                transcript.name,
                variant)
            return None

        sequence_key = ReferenceSequenceKey.from_variant_and_transcript(
            variant=variant,
            transcript=transcript,
            context_size=context_size)

        if sequence_key is None:
            logger.info(
                "No sequence key for variant %s on transcript %s",
                variant,
                transcript.name)
            return None

        return cls.from_variant_and_transcript_and_sequence_key(
            variant=variant,
            transcript=transcript,
            sequence_key=sequence_key)


def reading_frame_to_offset(reading_frame_at_start_of_sequence):
    """
    Given a reading frame (how many nucleotides into a codon) at the
    start of a cDNA sequence, return the number of nucleotides which need
    to be trimmed to start on a complete codon.

    Parameters
    ----------

    reading_frame_at_start_of_sequence : int

    Returns an int
    """
    if reading_frame_at_start_of_sequence < 0:
        raise ValueError("Reading frame can't be negative: %d" % (
            reading_frame_at_start_of_sequence,))
    elif reading_frame_at_start_of_sequence > 2:
        raise ValueError("Reading frame must be within 0 and 2, not %d" % (
            reading_frame_at_start_of_sequence,))
    # If we're 1 nucleotide into the codon then we need to shift
    # over two more to restore the ORF. Likewise, if we're 2 nucleotides in
    # then we have to shift over one more.
    return (3 - reading_frame_at_start_of_sequence) % 3
