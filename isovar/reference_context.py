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
from collections import namedtuple, OrderedDict, defaultdict
import logging

from skbio import DNA

from .effect_prediction import reference_transcripts_for_variant
from .variant_helpers import interbase_range_affected_by_variant_on_transcript
from .dataframe_builder import DataFrameBuilder


logger = logging.getLogger(__name__)


##########################
#
# SequenceKey
# -----------
#
# Used to identify and group the distinct sequences occurring  on a set of
# transcripts overlapping a variant locus
#
##########################

SequenceKey = namedtuple(
    "SequenceKey", [
        "strand",
        "sequence_before_variant_locus",
        "sequence_at_variant_locus",
        "sequence_after_variant_locus"
    ]
)


##########################
#
# SequenceKeyWithReadingFrame
# ---------------------------
#
# Includes all the fields of a SequenceKey, but also includes a reading frame
# at the start of the reference sequence.
#
#
##########################

SequenceKeyWithReadingFrame = namedtuple(
    "ReferenceContext", SequenceKey._fields + (
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
        "amino_acids_before_variant",
    )
)


##########################
#
# ReferenceContext
# ----------------
#
# Includes all the fields of SequenceKeyWithReadingFrame in addition to which
# variant we're examining and all transcripts overlapping that variant
# which produced this particular sequence context and reading frame.
#
##########################

ReferenceContext = namedtuple(
    "ReferenceContext",
    SequenceKeyWithReadingFrame._fields + (
        "variant",
        "transcripts")
)

def variant_matches_reference_sequence(variant, ref_seq_on_transcript, strand):
    """
    Make sure that reference nucleotides we expect to see on the reference
    transcript from a variant are the same ones we encounter.
    """
    if strand == "-":
        ref_seq_on_transcript = str(
            DNA(ref_seq_on_transcript).reverse_complement())

    return ref_seq_on_transcript == variant.ref

def sequence_key_for_variant_on_transcript(variant, transcript, context_size):
    """
    Extracts the reference sequence around a variant locus on a particular
    transcript.

    Parameters
    ----------
    variant : varcode.Variant

    transcript : pyensembl.Transcript

    context_size : int

    Returns SequenceKey object with the following fields:
        - strand
        - sequence_before_variant_locus
        - sequence_at_variant_locus
        - sequence_after_variant_locus

    Can also return None if Transcript lacks sufficiently long sequence
    """

    full_sequence = transcript.sequence

    if full_sequence is None:
        logger.warn(
            "Expected transcript %s (overlapping %s) to have sequence",
                transcript,
                variant)
        return None

    if len(full_sequence) < 6:
        # need at least 6 nucleotides for a start and stop codon
        logger.warn(
            "Sequence of %s (overlapping %s) too short: %d",
                transcript,
                variant,
                len(full_sequence))
        return None

    # get the interbase range of offsets which capture all reference
    # bases modified by the variant
    variant_start_offset, variant_end_offset = \
        interbase_range_affected_by_variant_on_transcript(
            variant=variant,
            transcript=transcript)

    logger.info("Interbase offset range on %s for variant %s = %d:%d",
        transcript,
        variant,
        variant_start_offset,
        variant_end_offset)

    prefix = full_sequence[
        max(0, variant_start_offset - context_size):
        variant_start_offset]

    suffix = full_sequence[
        variant_end_offset:
        variant_end_offset + context_size]

    ref_nucleotides_at_variant = full_sequence[
        variant_start_offset:variant_end_offset]
    if not variant_matches_reference_sequence(
            variant=variant,
            strand=transcript.strand,
            ref_seq_on_transcript=ref_nucleotides_at_variant):
        # TODO: once we're more confident about other logic in isovar,
        # change this to a warning and return None to allow for modest
        # differences between reference sequence patches, since
        # GRCh38.p1 may differ at some positions from GRCh38.p5
        raise ValueError(
            "Wrong reference sequence for variant %s on transcript %s" % (
                variant,
                transcript))
    return SequenceKey(
        strand=transcript.strand,
        sequence_before_variant_locus=prefix,
        sequence_at_variant_locus=ref_nucleotides_at_variant,
        sequence_after_variant_locus=suffix)


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


def sequence_key_with_reading_frame_for_variant_on_transcript(
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
                transcript,
                variant)
        return None

    if not transcript.contains_stop_codon:
        logger.info(
            "Expected transcript %s for variant %s to have stop codon",
                transcript,
                variant)
        return None

    if not transcript.protein_sequence:
        logger.info(
            "Expected transript %s for variant %s to have protein sequence",
                transcript,
                variant)
        return None

    sequence_key = sequence_key_for_variant_on_transcript(
        variant=variant,
        transcript=transcript,
        context_size=context_size)

    if sequence_key is None:
        logger.info("No sequence key for variant %s on transcript %s",
            variant,
            transcript)
        return None

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
                transcript,
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
    return SequenceKeyWithReadingFrame(
        strand=sequence_key.strand,
        sequence_before_variant_locus=sequence_key.sequence_before_variant_locus,
        sequence_at_variant_locus=sequence_key.sequence_at_variant_locus,
        sequence_after_variant_locus=sequence_key.sequence_after_variant_locus,
        offset_to_first_complete_codon=offset_to_first_complete_codon,
        contains_start_codon=contains_start_codon,
        overlaps_start_codon=overlaps_start_codon,
        contains_five_prime_utr=contains_five_prime_utr,
        amino_acids_before_variant=amino_acids_before_variant)

def sort_key_decreasing_max_length_transcript_cds(reference_context):
    """
    Used to sort a sequence of ReferenceContext objects by the longest CDS
    in each context's list of transcripts.
    """
    return -max(len(t.coding_sequence) for t in reference_context.transcripts)

def reference_contexts_for_variant(
        variant,
        context_size,
        transcript_id_whitelist=None):
    """
    variant : varcode.Variant

    context_size : int
        Max of nucleotides to include to the left and right of the variant
        in the context sequence.

    transcript_id_whitelist : set, optional
        If given, then only consider transcripts whose IDs are in this set.

    Returns list of ReferenceContext objects, sorted by maximum length of
    coding sequence of any supporting transcripts.
    """
    overlapping_transcripts = reference_transcripts_for_variant(
        variant=variant,
        transcript_id_whitelist=transcript_id_whitelist)

    # dictionary mapping SequenceKeyWithReadingFrame keys to list of
    # transcript objects
    sequence_groups = defaultdict(list)

    for transcript in overlapping_transcripts:
        sequence_key_with_reading_frame = \
            sequence_key_with_reading_frame_for_variant_on_transcript(
                variant=variant,
                transcript=transcript,
                context_size=context_size)
        if sequence_key_with_reading_frame is not None:
            sequence_groups[sequence_key_with_reading_frame].append(transcript)

    reference_contexts = [
        ReferenceContext(
            strand=key.strand,
            sequence_before_variant_locus=key.sequence_before_variant_locus,
            sequence_at_variant_locus=key.sequence_at_variant_locus,
            sequence_after_variant_locus=key.sequence_after_variant_locus,
            offset_to_first_complete_codon=key.offset_to_first_complete_codon,
            contains_start_codon=key.contains_start_codon,
            overlaps_start_codon=key.overlaps_start_codon,
            contains_five_prime_utr=key.contains_five_prime_utr,
            amino_acids_before_variant=key.amino_acids_before_variant,
            variant=variant,
            transcripts=matching_transcripts)
        for (key, matching_transcripts) in sequence_groups.items()
    ]
    reference_contexts.sort(key=sort_key_decreasing_max_length_transcript_cds)
    return reference_contexts

def reference_contexts_for_variants(
        variants,
        context_size,
        transcript_id_whitelist=None):
    """
    Extract a set of reference contexts for each variant in the collection.

    Parameters
    ----------
    variants : varcode.VariantCollection

    context_size : int
        Max of nucleotides to include to the left and right of the variant
        in the context sequence.

    transcript_id_whitelist : set, optional
        If given, then only consider transcripts whose IDs are in this set.

    Returns a dictionary from variants to lists of ReferenceContext objects,
    sorted by max coding sequence length of any transcript.
    """
    result = OrderedDict()
    for variant in variants:
        result[variant] = reference_contexts_for_variant(
            variant=variant,
            context_size=context_size,
            transcript_id_whitelist=transcript_id_whitelist)
    return result

def variants_to_reference_contexts_dataframe(
        variants,
        context_size,
        transcript_id_whitelist=None):
    """
    Given a collection of variants, find all reference sequence contexts
    around each variant.

    Parameters
    ----------
    variants : varcode.VariantCollection

    context_size : int
        Max of nucleotides to include to the left and right of the variant
        in the context sequence.

    transcript_id_whitelist : set, optional
        If given, then only consider transcripts whose IDs are in this set.

    Returns a DataFrame with {"chr", "pos", "ref", "alt"} columns for variants,
    as well as all the fields of ReferenceContext.
    """

    df_builder = DataFrameBuilder(
        ReferenceContext,
        exclude=["variant"],
        converters=dict(transcripts=lambda ts: ";".join(t.name for t in ts)),
        extra_column_fns={
            "gene": lambda variant, _: ";".join(variant.gene_names),
        })
    for variant, reference_contexts in reference_contexts_for_variants(
            variants=variants,
            context_size=context_size,
            transcript_id_whitelist=transcript_id_whitelist).items():
        df_builder.add_many(variant, reference_contexts)
    return df_builder.to_dataframe()
