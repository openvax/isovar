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

from .logging import create_logger
from .effect_prediction import reference_transcripts_for_variant

logger = create_logger(__name__)


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
# ReferenceContext
# ----------------
#
# Includes all the fields of a SequenceKey, in addition to which variant we're
# examining, all transcripts overlapping that variant which matched this
# particular sequence context
#
##########################

ReferenceContext = namedtuple(
    "ReferenceContext", SequenceKey._fields + ("variant", "transcripts")
)


##########################
#
# ReferenceContextWithORF
# ----------------
#
# Same fields as ReferenceContext but also includes a reading frame at the
# start of the reference sequence.
#
##########################

ReferenceContextWithORF = namedtuple(
    "ReferenceContextWithORF", ReferenceContext._fields + (
        "reading_frame_start_of_sequence",
        "first_codon_offset"
    )
)


def sequence_key_for_variant_on_transcript(variant, transcript, context_size):
    """
    Extracts the reference sequence around a variant locus on a particular
    transcript.

    Parameters
    ----------
    variant : varcode.Variant

    transcript : pyensembl.Transcript

    context_size : int

    Returns tuple with the following fields:
        - strand
        - sequence before variant locus
        - sequence at variant locus
        - sequence after variant locus

    Can also return None if Transcript lacks start codon or sequence
    """

    full_sequence = transcript.sequence

    if full_sequence is None:
        logger.warn(
            "Expected transcript %s (overlapping %s) to have sequence" % (
                transcript,
                variant))
        return None

    if len(full_sequence) < 6:
        # need at least 6 nucleotides for a start and stop codon
        logger.warn(
            "Sequence of %s (overlapping %s) too short: %d" % (
                transcript,
                variant,
                len(full_sequence)))
        return None

    if not transcript.contains_start_codon:
        logger.warn(
            "Expected transcript %s (overlapping %s)to have start codon" % (
                transcript,
                variant))
        return None

    start_codon_idx = min(transcript.start_codon_spliced_offsets)

    variant_on_transcript_start, variant_on_transcript_end = \
        interbase_range_affected_by_variant_on_transcript(
            variant=variant,
            transcript=transcript)

    if variant.is_insertion:
        # insertions don't actually affect the base referred to
        # by the start position of the variant, but rather the
        # variant gets inserted *after* that position
        context_sequence_end_idx = variant_on_transcript_start + 1
    else:
        # if not an insertion then the start offset of the variant
        # should actually point to a modified reference nucleotide,
        # so we should keep it as the upper bound of the "query" region
        # of the variant sequence
        context_sequence_end_idx = variant_on_transcript_start

    context_sequence_start_idx = max(
        start_codon_idx, context_sequence_end_idx - context_size)

    context_sequence_before_variant = full_sequences[
        context_sequence_start_idx:context_sequence_end_idx]

    n_ref_nucleotides_mutated = len(variant.ref)


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

    Returns set of ReferenceContext objects.
    """
    transcripts = reference_transcripts_for_variant(
        variant=variant,
        transcript_id_whitelist=transcript_id_whitelist)

    sequence_groups = defaultdict(list)

    for transcript in transcripts:
        sequence_key = sequence_around_variant_on_transcript(
            variant=variant,
            transcript=transcript,
            context_size=context_size)

def reference_contexts_for_variants(
        variants,
        context_size,
        transcript_id_whitelist=None):
    """
    variants : varcode.VariantCollection

    context_size : int
        Max of nucleotides to include to the left and right of the variant
        in the context sequence.

    transcript_id_whitelist : set, optional
        If given, then only consider transcripts whose IDs are in this set.

    Returns a dictionary from variants to sets of ReferenceContext objects.
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
    pass

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

def split_reference_context_by_reading_frames(reference_context):
    """
    Converts a single ReferenceContext into potentially multiple
    ReferenceContextWithReadingFrame objects.

    Parameters
    ----------
    reference_context : ReferenceContext

    A ReferenceContext conatains a reference cDNA sequence around a variant
    locus which is possibly shared by multiple transcript. However, there is a
    small possibility that the transcripts don't agree on the reading frame
    at the start of the sequence.

    To make sure we don't run into any strange edge-cases, we use this function
    to split a single ReferenceContext into potentially multiple
    ReferenceContext objects with disjoint sets of transcripts.
    """
    pass
