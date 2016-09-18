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


"""
Since multiple variant sequences can translate to the same amino acid sequence,
this module aggregates equivalent Translation objects into a single
ProteinSequence.
"""

from __future__ import print_function, division, absolute_import
from collections import namedtuple, defaultdict

from .default_parameters import (
    MIN_TRANSCRIPT_PREFIX_LENGTH,
    MAX_REFERENCE_TRANSCRIPT_MISMATCHES,
    PROTEIN_SEQUENCE_LENGTH,
    MAX_PROTEIN_SEQUENCES_PER_VARIANT,
    MIN_READS_SUPPORTING_VARIANT_CDNA_SEQUENCE,
)
from .dataframe_builder import dataframe_from_generator
from .translation import translate_variant_reads, TranslationKey
from .read_helpers import group_reads_by_allele
from .variant_helpers import trim_variant

##########################
#
# ProteinSequence
# ---------------
#
# Translated amino acid sequence aggregated across possibly multiple
# VariantSequence and ReferenceContext objects (e.g. imagine two distinct
# sequences which contain synonymous codons).
#
# This is the final result of the isovar variant->expressed peptide pipeline.
#
##########################

ProteinSequence = namedtuple(
    "ProteinSequence",
    TranslationKey._fields + (
        # list of all the Translation objects which support this distinct
        # amino acid sequence
        "translations",
        # number of reads overlapping the variant locus supporting any allele,
        # including the reference, alt, or anything else
        "overlapping_reads",
        # number of reads overlapping this locus which support the reference
        # allele
        "ref_reads",
        # total number of reads at the locus which contained the variant
        # nucleotides, even if they supported other phased sequences
        "alt_reads",
        # number of unique read names from all the VariantSequence objects
        # from each translation
        "alt_reads_supporting_protein_sequence",
        # how many reference transcripts overlap the variant locus?
        "transcripts_overlapping_variant",
        # how many reference transcripts were used to establish the
        # reading frame for this protein sequence
        "transcripts_supporting_protein_sequence",
        # name of gene of the reference transcripts used in Translation
        # objects
        "gene",
    ))


def to_translation_key(x):
    """
    If a namedtuple has a superset of the fields of a TranslationKey,
    project it down to create a TranslationKey object.
    """
    values_dict = {
        name: getattr(x, name)
        for name in TranslationKey._fields
    }
    return TranslationKey(**values_dict)


def group_translations(translations):
    """
    Returns a dictionary mapping TranslationKey object to a list of
    Translations with the same protein sequence
    """
    equivalent_translations_dict = defaultdict(list)
    for translation in translations:
        equivalent_translations_dict[to_translation_key(translation)].append(translation)
    return equivalent_translations_dict

def translation_key_to_protein_sequence(translation_key, **extra_kwargs):
    """
    Create a ProteinSequence object from a TranslationKey, and extra fields
    which must be supplied in extra_kwargs.
    """
    values_dict = {}
    for name in TranslationKey._fields:
        values_dict[name] = getattr(translation_key, name)
    values_dict.update(extra_kwargs)
    return ProteinSequence(**values_dict)

def summarize_translations(translations):
    """
    Summarize a collection of Translation objects into three values:
        1) List of unique reads supporting underlying variant sequences
        2) Set of unique transcript names for establishing reading frames of the
           translations.
        3) Set of unique gene names for all transcripts used by translations.
    """
    read_name_to_reads = {}
    gene_names = set([])
    transcript_ids = set([])
    for translation in translations:
        for read in translation.variant_sequence.reads:
            read_name_to_reads[read.name] = read
        for transcript in translation.reference_context.transcripts:
            transcript_ids.add(transcript.id)
            gene_names.add(transcript.gene.name)
    unique_reads = list(read_name_to_reads.values())
    return unique_reads, transcript_ids, gene_names

def protein_sequence_sort_key(protein_sequence):
    """
    Sort protein sequences lexicographically by three criteria:
        - number of unique supporting reads
        - minimum mismatch versus a supporting reference transcript
        - number of supporting reference transcripts
    """
    return (
        len(protein_sequence.alt_reads_supporting_protein_sequence),
        min(
            t.number_mismatches
            for t in protein_sequence.translations),
        len(protein_sequence.transcripts_supporting_protein_sequence)
    )

def sort_protein_sequences(protein_sequences):
    """
    Sort protein sequences in decreasing order of priority
    """
    return list(
        sorted(protein_sequences, key=protein_sequence_sort_key, reverse=True))

def reads_generator_to_protein_sequences_generator(
        variant_and_overlapping_reads_generator,
        transcript_id_whitelist=None,
        protein_sequence_length=PROTEIN_SEQUENCE_LENGTH,
        min_reads_supporting_cdna_sequence=MIN_READS_SUPPORTING_VARIANT_CDNA_SEQUENCE,
        min_transcript_prefix_length=MIN_TRANSCRIPT_PREFIX_LENGTH,
        max_transcript_mismatches=MAX_REFERENCE_TRANSCRIPT_MISMATCHES,
        max_protein_sequences_per_variant=MAX_PROTEIN_SEQUENCES_PER_VARIANT):
    """"
    Translates each coding variant in a collection to one or more
    Translation objects, which are then aggregated into equivalent
    ProteinSequence objects.

    Parameters
    ----------
    transcript_id_whitelist : set, optional
        If given, expected to be a set of transcript IDs which we should use
        for determining the reading frame around a variant. If omitted, then
        try to use all overlapping reference transcripts.

    protein_sequence_length : int
        Try to translate protein sequences of this length, though sometimes
        we'll have to return something shorter (depending on the RNAseq data,
        and presence of stop codons).

    min_reads_supporting_cdna_sequence : int
        Drop variant sequences supported by fewer than this number of reads.

    min_transcript_prefix_length : int
        Minimum number of bases we need to try matching between the reference
        context and variant sequence.

    max_transcript_mismatches : int
        Don't try to determine the reading frame for a transcript if more
        than this number of bases differ.

    max_protein_sequences_per_variant : int
        Number of protein sequences to return for each ProteinSequence


    Yields pairs of a Variant and a list of ProteinSequence objects
    """

    for (variant, overlapping_reads) in variant_and_overlapping_reads_generator:
        overlapping_transcript_ids = [
            t.id
            for t in variant.transcripts
            if t.is_protein_coding
        ]
        _, ref, alt = trim_variant(variant)
        overlapping_reads = list(overlapping_reads)
        reads_grouped_by_allele = group_reads_by_allele(overlapping_reads)

        ref_reads = reads_grouped_by_allele.get(ref, [])
        alt_reads = reads_grouped_by_allele.get(alt, [])

        translations = translate_variant_reads(
            variant=variant,
            variant_reads=alt_reads,
            transcript_id_whitelist=transcript_id_whitelist,
            protein_sequence_length=protein_sequence_length,
            min_reads_supporting_cdna_sequence=min_reads_supporting_cdna_sequence,
            min_transcript_prefix_length=min_transcript_prefix_length,
            max_transcript_mismatches=max_transcript_mismatches)

        protein_sequences = []
        for (key, equivalent_translations) in group_translations(translations).items():
            # get the variant read names, transcript IDs and gene names for
            # protein sequence we're about to construct
            alt_reads_supporting_protein_sequence, group_transcript_ids, group_gene_names = \
                summarize_translations(equivalent_translations)

            protein_sequence = translation_key_to_protein_sequence(
                translation_key=key,
                translations=equivalent_translations,
                overlapping_reads=overlapping_reads,
                alt_reads=alt_reads,
                ref_reads=ref_reads,
                alt_reads_supporting_protein_sequence=alt_reads_supporting_protein_sequence,
                transcripts_supporting_protein_sequence=group_transcript_ids,
                transcripts_overlapping_variant=overlapping_transcript_ids,
                gene=list(group_gene_names))
            protein_sequences.append(protein_sequence)

        # sort protein sequences before returning the top results
        protein_sequences = sort_protein_sequences(protein_sequences)

        yield variant, protein_sequences[:max_protein_sequences_per_variant]

def protein_sequences_generator_to_dataframe(variant_and_protein_sequences_generator):
    """
    Given a generator which yields (Variant, [ProteinSequence]) elements,
    returns a pandas.DataFrame
    """
    return dataframe_from_generator(
        element_class=ProteinSequence,
        variant_and_elements_generator=variant_and_protein_sequences_generator,
        converters=dict(
            gene=lambda x: ";".join(x)))
