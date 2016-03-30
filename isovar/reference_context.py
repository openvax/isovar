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

import pandas as pd

from .logging import create_logger
from .effect_prediction import reference_transcripts_for_variant
from .variant_helpers import interbase_range_affected_by_variant_on_transcript

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
        "reading_frame",
        "offset_to_first_codon"
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
    SequenceKeyWithReadingFrame._fields + ("variant", "transcripts")
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

    # get the interbase range of offsets which capture all reference
    # bases modified by the variant
    variant_start_offset, variant_end_offset = \
        interbase_range_affected_by_variant_on_transcript(
            variant=variant,
            transcript=transcript)
    print(variant_start_offset, variant_end_offset)
    prefix = full_sequence[
        max(0, variant_start_offset - context_size):
        variant_start_offset]

    suffix = full_sequence[
        variant_end_offset:
        variant_end_offset + context_size]

    variant_nucleotides = full_sequence[variant_start_offset:variant_end_offset]

    return SequenceKey(
        strand=transcript.strand,
        sequence_before_variant_locus=prefix,
        sequence_at_variant_locus=variant_nucleotides,
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

    Returns SequenceKeyWithReadingFrame object with the following fields:
        - strand
        - sequence_before_variant_locus
        - sequence_at_variant_locus
        - sequence_after_variant_locus
        - reading_frame
        - offset_to_first_codon

    Can also return None if Transcript lacks sufficiently long sequence or
    annotated start/stop codons.
    """
    if not transcript.contains_start_codon:
        logger.info(
            "Expected transcript %s for variant %s to have start codon" % (
                transcript,
                variant))
        return None

    if not transcript.contains_stop_codon:
        logger.info(
            "Expected transcript %s for variant %s to have stop codon" % (
                transcript,
                variant))
        return None

    sequence_key = sequence_key_for_variant_on_transcript(
        variant=variant,
        transcript=transcript,
        context_size=context_size)

    if sequence_key is None:
        logger.info("No sequence key for variant %s on transcript %s" % (
            variant,
            transcript))
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
    if variant_start_offset <= start_codon_idx + 3:
        logger.info(
            "Skipping transcript %s for variant %s, must be after start codon" % (
                transcript,
                variant))
        return None

    stop_codon_idx = min(transcript.stop_codon_spliced_offsets)

    # skip variants which affect the 3' UTR of the transcript since
    # they don't have obvious coding effects on the protein sequence
    if variant_start_offset >= stop_codon_idx + 3:
        logger.info(
            "Skipping transcript %s for variant %s, occurs in 3' UTR" % (
                transcript,
                variant))
        return None

    n_prefix = len(sequence_key.sequence_before_variant_locus)
    prefix_start_idx = variant_start_offset - n_prefix

    reading_frame = (prefix_start_idx - start_codon_idx) % 3
    offset_to_first_codon = reading_frame_to_offset(reading_frame)
    return SequenceKeyWithReadingFrame(
        strand=sequence_key.strand,
        sequence_before_variant_locus=sequence_key.sequence_before_variant_locus,
        sequence_at_variant_locus=sequence_key.sequence_at_variant_locus,
        sequence_after_variant_locus=sequence_key.sequence_after_variant_locus,
        reading_frame=reading_frame,
        offset_to_first_codon=offset_to_first_codon)

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
            strand=sequence_key.strand,
            sequence_before_variant_locus=sequence_key.sequence_before_variant_locus,
            sequence_at_variant_locus=sequence_key.sequence_at_variant_locus,
            sequence_after_variant_locus=sequence_key.sequence_after_variant_locus,
            reading_frame=sequence_key.reading_frame,
            offset_to_first_codon=sequence_key.offset_to_first_codon,
            variant=variant,
            transcripts=matching_transcripts)
        for (sequence_key, matching_transcripts) in sequence_groups.items()
    ]
    reference_contexts.sort(key=sort_key_decreasing_max_length_transcript_cds)
    return reference_contexts

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
    columns = [
        ("chr", []),
        ("pos", []),
        ("ref", []),
        ("alt", []),
    ]
    for field in ReferenceContext._fields:
        columns.append((field, []))
    columns_dict = OrderedDict(columns)

    for variant, reference_context in reference_contexts_for_variants(
            variants=variants,
            context_size=context_size,
            transcript_id_whitelist=transcript_id_whitelist):
        columns_dict["chr"].append(variant.contig)
        columns_dict["pos"].append(variant.original_start)
        columns_dict["ref"].append(variant.original_ref)
        columns_dict["alt"].append(variant.original_alt)
        for field in ReferenceContext._fields:
            columns_dict[field].append(getattr(reference_context, field))

    return pd.DataFrame(columns_dict)
