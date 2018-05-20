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

from __future__ import print_function, division, absolute_import

from .value_object import ValueObject
from .dna import reverse_complement_dna
from .variant_helpers import interbase_range_affected_by_variant_on_transcript
from .logging import get_logger

logger = get_logger(__name__)


class ReferenceSequenceKey(ValueObject):
    """
    Used to identify and group the distinct sequences occurring on a set of
    transcripts overlapping a variant locus.
    """
    __slots__ = [
        "strand",
        "sequence_before_variant_locus",
        "sequence_at_variant_locus",
        "sequence_after_variant_locus",
    ]

    def __init__(
            self,
            strand,
            sequence_before_variant_locus,
            sequence_at_variant_locus,
            sequence_after_variant_locus):
        if strand not in {'+', '-'}:
            raise ValueError("Invalid strand: '%s'" % strand)
        self.strand = strand
        self.sequence_before_variant_locus = sequence_before_variant_locus
        self.sequence_at_variant_locus = sequence_at_variant_locus
        self.sequence_after_variant_locus = sequence_after_variant_locus

    @classmethod
    def from_variant_and_transcript(
            cls, variant, transcript, context_size):
        """
        Extracts the reference sequence around a variant locus on a particular
        transcript.

        Parameters
        ----------
        variant : varcode.Variant

        transcript : pyensembl.Transcript

        context_size : int

        Returns SequenceKey object.

        Can also return None if Transcript lacks sufficiently long sequence
        """

        full_transcript_sequence = transcript.sequence

        if full_transcript_sequence is None:
            logger.warn(
                "Expected transcript %s (overlapping %s) to have sequence",
                transcript.name,
                variant)
            return None

        # get the interbase range of offsets which capture all reference
        # bases modified by the variant
        variant_start_offset, variant_end_offset = \
            interbase_range_affected_by_variant_on_transcript(
                variant=variant,
                transcript=transcript)

        reference_cdna_at_variant = full_transcript_sequence[
            variant_start_offset:variant_end_offset]

        if not variant_matches_reference_sequence(
                variant=variant,
                strand=transcript.strand,
                ref_seq_on_transcript=reference_cdna_at_variant):
            logger.warn(
                "Variant %s doesn't match reference sequence on transcript %s: "
                "may span splice junction",
                variant,
                transcript)
            return None

        if len(full_transcript_sequence) < 6:
            # need at least 6 nucleotides for a start and stop codon
            logger.warn(
                "Sequence of %s (overlapping %s) too short: %d",
                transcript,
                variant,
                len(full_transcript_sequence))
            return None

        logger.info(
            "Interbase offset range on %s for variant %s = %d:%d",
            transcript.name,
            variant,
            variant_start_offset,
            variant_end_offset)

        reference_cdna_before_variant = full_transcript_sequence[
            max(0, variant_start_offset - context_size):
            variant_start_offset]

        reference_cdna_after_variant = full_transcript_sequence[
            variant_end_offset:
            variant_end_offset + context_size]

        return ReferenceSequenceKey(
            strand=transcript.strand,
            sequence_before_variant_locus=reference_cdna_before_variant,
            sequence_at_variant_locus=reference_cdna_at_variant,
            sequence_after_variant_locus=reference_cdna_after_variant)


def variant_matches_reference_sequence(variant, ref_seq_on_transcript, strand):
    """
    Make sure that reference nucleotides we expect to see on the reference
    transcript from a variant are the same ones we encounter.
    """
    if strand == "-":
        ref_seq_on_transcript = reverse_complement_dna(ref_seq_on_transcript)
    return ref_seq_on_transcript == variant.ref
