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
import logging

from skbio import DNA
from varcode import EffectCollection

def predicted_coding_effects_with_mutant_sequence(
        variant,
        transcript_id_whitelist=None):
    """
    For a given variant, return the set of predicted mutation effects
    on transcripts where this variant results in a predictable non-synonymous
    change to the protein sequence.

    Parameters
    ----------
    variant : varcode.Variant

    transcript_id_whitelist : set
        Filter effect predictions to only include these transcripts

    Returns a varcode.EffectCollection object
    """

    effects = []
    for transcript in variant.transcripts:
        if not transcript.complete:
            logging.info(
                "Skipping transcript %s for variant %s because it's incomplete" % (
                    transcript,
                    variant))
            continue

        if transcript_id_whitelist and transcript.id not in transcript_id_whitelist:
            logging.info(
                "Skipping transcript %s for variant %s because it's not one of %d allowed" % (
                    transcript,
                    variant,
                    len(transcript_id_whitelist)))
            continue
        effects.append(variant.effect_on_transcript(transcript))

    effects = EffectCollection(effects)

    n_total_effects = len(effects)
    logging.info("Predicted %d effects for variant %s" % (
        n_total_effects,
        variant))

    nonsynonymous_coding_effects = effects.drop_silent_and_noncoding()
    logging.info(
        "Keeping %d/%d non-synonymous coding effects for %s" % (
            len(nonsynonymous_coding_effects),
            n_total_effects,
            variant))

    usable_effects = [
        effect
        for effect in nonsynonymous_coding_effects
        if effect.mutant_protein_sequence is not None
    ]
    logging.info(
        "Keeping %d/%d effects with predictable AA sequences for %s" % (
            len(usable_effects),
            len(nonsynonymous_coding_effects),
            variant))
    return usable_effects

def reference_transcripts_for_variant(variant, transcript_id_whitelist=None):
    """
    For a given variant, find all the transcripts which overlap the
    variant and for which it has a predictable effect on the amino acid
    sequence of the protein.
    """
    predicted_effects = predicted_coding_effects_with_mutant_sequence(
        variant=variant,
        transcript_id_whitelist=transcript_id_whitelist)
    return [effect.transcript for effect in predicted_effects]


ReferenceContext = namedtuple(
    "ReferenceContext",
    [
        "strand",
        "transcript_ids",
        "transcript_names",
        "gene",
        "reference_cdna_sequence_around_variant",
        "reading_frame_at_start_of_context",
        "number_cdna_bases_before_variant",
        "number_cdna_bases_after_variant",
        "reference_protein_sequence_around_variant",
    ])

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
    # First, we must convert it to a stranded offset on the mRNA.
    # We must then adjust this offset into a start:end range depending
    # on which kind of variant we're dealing with and whether the
    # transcript is actually on the backward strand.
    forward_strand_offset = transcript.spliced_offset(variant.start)

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
                logging.info(
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

def variants_to_reference_contexts(
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

    Returns a dictionary from variants to ReferenceContext objects.
    """
    result = OrderedDict()
    for variant in variants:
        effects = predicted_coding_effects_with_mutant_sequence(
            variant=variant,
            transcript_id_whitelist=transcript_id_whitelist)
        transcripts = [effect.transcript for effect in effects]
        for transcript in transcripts:

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