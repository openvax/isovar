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

from skbio import DNA

from .common import create_logger

logger = create_logger(__name__)


MAX_TRANSCRIPT_MISMATCHES = 2
MIN_TRANSCRIPT_PREFIX_LENGTH = 10

# information related to the translation of a RNA sequence context
# in a reading frame determined by a particular reference transcript
TranslationFromReferenceORF = namedtuple(
    "TranslationFromReferenceORF",
    [
        "cdna_prefix",
        "cdna_variant",
        "cdna_suffix",
        "transcript_id",
        "transcript_name",
        "number_transcript_sequence_mismatches",
        "reading_frame_at_start_of_cdna_sequence",
        "transcript_sequence_before_variant",
        "variant_protein_sequence",
        "reference_protein_sequence",
        "fragment_aa_start_offset_in_protein",
        "fragment_aa_end_offset_in_protein"
        "variant_aa_start_offset_in_fragment",
        "variant_aa_end_offset_in_fragment",
    ])

class TranslationOrFailure(object):
    def __init__(self, translation_info=None, failure_message=None):
        self.translation_info = translation_info
        self.failure_message = failure_message

        none_count = 0
        none_count += (self.translation_info is None)
        none_count += (self.failure_message is None)
        if none_count == 0:
            raise ValueError(
                "One of translation_info or failure_message must be not None")
        elif none_count == 2:
            raise ValueError(
                "One of translation_info and failure_message must have a valid value")


def _translate(
        cdna_sequence_before_variant,
        cdna_sequence_of_variant,
        cdna_sequence_after_variant,
        transcript_id,
        transcript_reading_frame_at_start_of_context_sequence,
        transcript_full_cdna,
        transcript_full_protein,
        variant_base1_offset_in_transcript):
    """
    Try to  translate the cDNA sequence context around a variant using the
    reading frame of a transcript.

    We're assuming that the sequence of the reference transcript up to the
    variant largely matches the sequence we detected from RNA.
    """
    if transcript_reading_frame_at_start_of_context_sequence == 1:
        # if we're 1 nucleotide into the codon then we need to shift
        # over two more to restore the ORF
        orf_offset = 2
    elif transcript_reading_frame_at_start_of_context_sequence == 2:
        orf_offset = 1
    else:
        orf_offset = 0

    logger.info("ORF offset into sequence %s_%s_%s from transcript %s: %d" % (
        cdna_sequence_before_variant,
        cdna_sequence_of_variant,
        cdna_sequence_after_variant,
        transcript_id,
        orf_offset))

    # translate the variant cDNA sequence we detected from spanning reads
    # using the ORF offset from the current reference transcript
    combined_variant_cdna_sequence = DNA(
        cdna_sequence_before_variant + cdna_sequence_of_variant + cdna_sequence_after_variant)

    variant_protein_fragment_sequence = str(
        combined_variant_cdna_sequence[orf_offset:].translate())

    logger.info("Combined variant cDNA sequence %s translates to protein %s" % (
        combined_variant_cdna_sequence,
        variant_protein_fragment_sequence))

    # translate the reference sequence using the given ORF offset,
    # we can probably sanity check this by making sure it matches a
    # substring of the transcript.protein_sequence field
    combined_transcript_sequence = DNA(
        transcript_full_cdna[
            query_sequence_start_idx:
            query_sequence_start_idx + len(combined_variant_cdna_sequence)])

    transcript_protein_fragment_sequence = str(
        combined_transcript_sequence[orf_offset].translate())

    fragment_aa_start_offset_in_protein = (
        query_sequence_start_idx - start_codon_idx) // 3
    fragment_aa_end_offset_in_protein = (
        query_sequence_end_idx - start_codon_idx) // 3 + 1

    n_prefix_codons = n_prefix_nucleotides // 3

    if variant_is_deletion and variant_is_frameshift:
        # after a deletion that shifts the reading frame,
        # there is a partial codon left, which may be different from
        # the original codon at that location. Additionally,
        # the reading frame for all subsequent codons will be different
        # from the reference
        variant_codons_start_offset_in_fragment = None

    # TODO: change variant_codons_start_offset_in_fragment in the
    # TranslationFromReferenceORF to
    # - variant_amino_acids_start_offset_in_fragment
    # - variant_amino_acids_end_offset_in_fragment
    # which are computed from _codons_ offsets by checking which amino
    # acids differ from the reference sequence(s).
    # If a mutated sequence ends up having no differences then
    # we should skip it as a synonymous mutation.

    # TODO #2:
    # Instead of doing this once per sequence/transcript pair, we
    # should precompute the reference transcript sequences and the
    # ORFs they imply and pass in a list of ReferenceTranscriptContext
    # objects with the following field:
    #   sequence_before_variant_locus : cDNA
    #   reading_frame_at_start_of_sequence : int
    #   transcript_ids : str set
    #   transcript_names : str set
    #   gene : str

    variant_aa_start_offset_in_fragment = 0
    last_variant_codon = 0

    """
    # the number of non-mutated codons in the prefix (before the variant)
    # has to trim the ORF offset and then count up by multiples of 3


    aa_prefix = variant_protein_fragment_sequence[:n_prefix_codons]
    aa_variant = variant_protein_fragment_sequence[
        n_prefix_codons:n_prefix_codons + n_variant_codons]
    aa_suffix = variant_protein_fragment_sequence[
        n_prefix_codons + n_variant_codons:]

    assert aa_prefix + aa_variant + aa_suffix == variant_protein_fragment_sequence
                    aa_prefix=aa_prefix,
            aa_variant=aa_variant,
            aa_suffix=aa_suffix)
    """



def translate_variant_on_transcript(
        dna_sequence_prefix,
        dna_sequence_variant,
        dna_sequence_suffix,
        variant_is_insertion,
        variant_is_deletion,
        variant_is_frameshift,
        base1_variant_start_location,
        base1_variant_end_location,
        transcript,
        max_transcript_mismatches=MAX_TRANSCRIPT_MISMATCHES,
        min_transcript_prefix_length=MIN_TRANSCRIPT_PREFIX_LENGTH):
    if transcript.strand == "+":
        cdna_prefix = dna_sequence_prefix
        cdna_suffix = dna_sequence_suffix
        cdna_variant = dna_sequence_variant
        variant_in_transcript_idx = transcript.spliced_offset(
            base1_variant_start_location)
    else:
        # if the transcript is on the reverse strand then we have to
        # take the sequence PREFIX|VARIANT|SUFFIX
        # and take the complement of XIFFUS|TNAIRAV|XIFERP
        cdna_prefix = str(DNA(dna_sequence_suffix).reverse_complement())
        cdna_suffix = str(DNA(dna_sequence_prefix).reverse_complement())
        cdna_variant = str(
            DNA(dna_sequence_variant).reverse_complement())
        variant_in_transcript_idx = transcript.spliced_offset(
            base1_variant_end_location)

    if variant_is_insertion:
        # insertions don't actually affect the base referred to
        # by the start position of the variant, but rather the
        # variant gets inserted *after* that position
        query_sequence_end_idx = variant_in_transcript_idx + 1
    else:
        query_sequence_end_idx = variant_in_transcript_idx

    start_codon_idx = min(transcript.start_codon_spliced_offsets)

    if variant_in_transcript_idx < start_codon_idx + 3:
        logger.info(
            "Skipping %s because variant appears in 5' UTR" % (
                transcript))
        return None

    query_sequence_start_idx = query_sequence_end_idx - len(cdna_prefix)

    if query_sequence_start_idx < 0:
        logger.warn("Transcript %s not long enough for observed sequence" % (
            transcript))
        return None

    transcript_sequence_before_variant = transcript.sequence[
        query_sequence_start_idx:query_sequence_end_idx]

    assert len(transcript_sequence_before_variant) == len(cdna_prefix)

    if len(transcript_sequence_before_variant) < min_transcript_prefix_length:
        logger.info(
            "Skipping transcript %s because it does not have %d nucleotides before variant" % (
                transcript,
                min_transcript_prefix_length))
        return None

    n_mismatch_before_variant = sum(
        xi != yi
        for (xi, yi) in zip(
            transcript_sequence_before_variant, cdna_prefix))

    if n_mismatch_before_variant > max_transcript_mismatches:
        logger.info(
            "Skipping transcript %s, too many mismatching bases (%d)",
            transcript,
            n_mismatch_before_variant)
        return None

    transcript_reading_frame = (query_sequence_start_idx - start_codon_idx) % 3

    return TranslationFromReferenceORF(
        cdna_prefix=cdna_prefix,
        cdna_variant=cdna_variant,
        cdna_suffix=cdna_suffix,
        transcript_id=transcript.id,
        transcript_name=transcript.name,
        number_transcript_sequence_mismatches=n_mismatch,
        reading_frame_at_start_of_cdna_sequence=reading_frame,
        transcript_sequence_before_variant=transcript_sequence_before_variant,
        variant_protein_sequence=variant_protein_fragment_sequence,
        reference_protein_sequence=transcript_protein_fragment_sequence,
        fragment_aa_start_offset_in_protein=fragment_aa_start_offset_in_protein,
        fragment_aa_end_offset_in_protein=fragment_aa_end_offset_in_protein,
        variant_aa_start_offset_in_fragment=variant_aa_start_offset_in_fragment,
        variant_aa_end_offset_in_fragment=variant_aa_end_offset_in_fragment)


def translate_compatible_reading_frames(
        dna_sequence_prefix,
        dna_sequence_variant,
        dna_sequence_suffix,
        variant_is_insertion,
        variant_is_deletion,
        variant_is_frameshift,
        base1_variant_start_location,
        base1_variant_end_location,
        transcripts,
        max_transcript_mismatches=MAX_TRANSCRIPT_MISMATCHES,
        min_transcript_prefix_length=MIN_TRANSCRIPT_PREFIX_LENGTH):
    """
    Use the given reference transcripts to attempt to establish the ORF
    of a sequence context extract from RNA reads. It's expected that the
    sequence has been aligned to the reference genome potentially using its
    reverse-complement, and thus the sequence we are given is always from the
    positive strand of DNA.

    Parameters
    ----------
    dna_sequence_prefix : str
        Nucleotides before the variant.

    dna_sequence_variant : str
        Mutated nucleotides, should be empty for a deletion.

    dna_sequence_suffix : str
        Nucleotides after the variant

    variant_is_insertion : bool

    base1_variant_start_location : int
        For deletions and substitutions, this position is the first modified
        nucleotide. For insertions, this is the position before any inserted
        nucleotides.

    base1_variant_end_location : int
        For deletions and substitutions, this is the positions of the last
        affected reference nucleotide. For insertions, this is the location of
        the base after the insertion.

    transcripts : list of pyensembl.Transcript
        List of candidate reference transcripts from which we try to determine
        the ORF of a variant RNA sequence.

    max_transcript_mismatches : int
        Ignore transcripts with more than this number of mismatching nucleotides

    min_transcript_prefix_length : int
        Don't consider transcripts with less than this number of nucleotides
        before the variant position (setting this value to 0 will enable use
        of transcripts without 5' UTR)

    Returns list of TranslationFromReferenceORF objects.
    """
    return [

    ]
