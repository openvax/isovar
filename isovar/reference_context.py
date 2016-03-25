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


reference_context_fields = [
    "variant",
    "strand",
    "reference_cdna_before_variant",
    "refernece_cdna_at_variant",
    "reference_cdna_after_variant",
    "reference_protein_sequence_around_variant",
    "transcripts",
]

ReferenceContext = namedtuple(
    "ReferenceContext",
    reference_context_fields)

ReferenceContextWithReadingFrame = namedtuple(
    "ReferenceContextWithReadingFrame",
    reference_context_fields + [
        "reading_frame_start_of_sequence",
        "first_codon_offset"])

def interbase_range_affected_by_variant_on_transcript(variant, transcript):
    """
    Convert from a variant's position in global genomic coordinates on the
    forward strand to an interval of interbase offsets on a particular
    transcript's mRNA.

    Parameters
    ----------
    variant : varcode.Variant

    transcript : pyensembl.Transcript

    Assumes that the transcript overlaps the variant.

    Returns (start, end) tuple of offsets into the transcript's cDNA sequence
    which indicates which bases in the reference sequence are affected by a
    variant.

    Example:
        The insertion of "TTT" into the middle of an exon would result in an
        offset pair such as (100,100) since no reference bases are changed
        or deleted by an insertion.

        On the other hand, deletion the preceding "CGG" at that same locus could
        result in an offset pair such as (97, 100)
    """
    if variant.is_insertion:
        if transcript.strand == "+":
            # base-1 position of an insertion is the genomic nucleotide
            # before any inserted mutant nucleotides
            start_offset = [transcript.spliced_offset(variant.start)]
        else:
            # assuming that this transcript was considered to overlap
            # with the variant since the insertion happens inside
            # one of its exons (rather than simply immediately before
            # of after)
            start_offset = [transcript.spliced_offset(variant.start + 1)]
        # an insertion happens *between* two reference bases
        # so the start:end offsets coincide
        end_offset = start_offset
    else:
        # reference bases affected by substitution or deletion defined by
        # range starting at first affected base
        offsets = []
        assert len(variant.ref) > 0
        for dna_pos in range(variant.start, variant.start + len(variant.ref)):
            try:
                offsets.append(transcript.spliced_offset(dna_pos))
            except ValueError:
                logger.info(
                    "Couldn't find position %d from %s on exons of %s" % (
                        dna_pos,
                        variant,
                        transcript))
        if len(offsets) == 0:
            raise ValueError(
                "Couldn't find any exonic reference bases affected by %s on %s" % (
                    variant,
                    transcript))
        start_offset = min(offsets)
        end_offset = max(offsets) + 1
    return (start_offset, end_offset)

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

    sequence_groups = {}
    for transcript in transcripts:
        variant_on_transcript_start, variant_on_transcript_end = \
            interbase_range_affected_by_variant_on_transcript(
                variant=variant,
                transcript=transcript)
        start_codon_idx = min(transcript.start_codon_spliced_offsets)

        if variant.is_insertion:
            # insertions don't actually affect the base referred to
            # by the start position of the variant, but rather the
            # variant gets inserted *after* that position
            sequence_end_idx = variant_on_transcript_start + 1
        else:
            # if not an insertion then the start offset of the variant
            # should actually point to a modified reference nucleotide,
            # so we should keep it as the upper bound of the "query" region
            # of the variant sequence
            sequence_end_idx = variant_on_transcript_start

        sequence_start_idx = max(
            start_codon_idx, sequence_end_idx - context_size)

        sequence = transcript.sequence[sequence_start_idx:sequence_end_idx]


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
    return 3 - reading_frame_at_start_of_sequence

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
