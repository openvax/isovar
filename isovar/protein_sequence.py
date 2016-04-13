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
    PROTEIN_SEQUENCE_LEGNTH,
    MAX_PROTEIN_SEQUENCES_PER_VARIANT,
    MIN_READS_SUPPORTING_VARIANT_CDNA_SEQUENCE,
)
from .dataframe_builder import DataFrameBuilder
from .translation import translate_variants, TranslationKey

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
        "translations_count",
        # number of unique read names from all the VariantSequence objects
        # from each translation
        "supporting_variant_reads_count",
        # total number of reads at the locus which contained the variant
        # nucleotides, even if they supported other phased sequences
        "total_variant_reads",
        # how many reference transcripts were used to establish the
        # reading frame for this protein sequence
        "supporting_transcripts_count",
        # how many reference transcripts overlap the variant locus?
        "total_transcripts",
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
        1) Set of unique read names supporting underlying variant sequences
        2) Set of unique transcript names for establishing reading frames of the
           translations.
        3) Set of unique gene names for all transcripts used by translations.
    """
    read_names = set([])
    gene_names = set([])
    transcript_ids = set([])
    for translation in translations:
        for read_name in translation.variant_sequence.read_names:
            read_names.add(read_name)
        for transcript in translation.reference_context.transcripts:
            transcript_ids.add(transcript.id)
            gene_names.add(transcript.gene.name)
    return read_names, transcript_ids, gene_names

def variants_to_protein_sequences(
        variants,
        samfile,
        transcript_id_whitelist=None,
        protein_sequence_length=PROTEIN_SEQUENCE_LEGNTH,
        min_reads_supporting_rna_sequence=MIN_READS_SUPPORTING_VARIANT_CDNA_SEQUENCE,
        min_transcript_prefix_length=MIN_TRANSCRIPT_PREFIX_LENGTH,
        max_transcript_mismatches=MAX_REFERENCE_TRANSCRIPT_MISMATCHES,
        max_protein_sequences_per_variant=MAX_PROTEIN_SEQUENCES_PER_VARIANT):

    for (variant, translations) in translate_variants(
            variants=variants,
            samfile=samfile,
            transcript_id_whitelist=transcript_id_whitelist,
            protein_sequence_length=protein_sequence_length,
            min_reads_supporting_rna_sequence=min_reads_supporting_rna_sequence,
            min_transcript_prefix_length=min_transcript_prefix_length,
            max_transcript_mismatches=max_transcript_mismatches):

        # convert to list so we can traverse it twice, once to get TranslationKey
        # mappings for each Translation and a second time to collect all the
        # unique read names and transript IDs.
        translations = list(translations)

        # dictionary mapping TranslationKey object to a list of Translations
        # with the same protein sequence
        equivalent_translations_dict = defaultdict(list)

        for translation in translations:
            key = to_translation_key(translation)
            equivalent_translations_dict[key].append(translation)

        all_read_names, all_transcript_ids, _ = summarize_translations(translations)

        protein_sequences = []
        for (key, equivalent_translations) in equivalent_translations_dict.items():
            # get the variant read names, transcript IDs and gene names for
            # protein sequence we're about to construct
            group_read_names, group_transcript_ids, group_gene_names = \
                summarize_translations(equivalent_translations)

            protein_sequence = translation_key_to_protein_sequence(
                translation_key=key,
                translations=equivalent_translations,
                translations_count=len(equivalent_translations),
                supporting_variant_reads_count=len(group_read_names),
                total_variant_reads=len(all_read_names),
                supporting_transcripts_count=len(group_transcript_ids),
                total_transcripts=len(all_transcript_ids),
                gene=list(group_gene_names))
            protein_sequences.append(protein_sequence)
        yield variant, protein_sequences[:max_protein_sequences_per_variant]

def variants_to_protein_sequences_dataframe(*args, **kwargs):
    """
    Given a collection of variants and a SAM/BAM file of overlapping reads,
    returns a DataFrame with a row for each protein sequence.

    Takes the same parameters as variants_to_protein_sequences.
    """
    df_builder = DataFrameBuilder(
        ProteinSequence,
        exclude_field_names=["translations"])
    for (variant, protein_sequences) in variants_to_protein_sequences(*args, **kwargs):
        for protein_sequence in protein_sequences:
            df_builder.add(variant, protein_sequence)
    return df_builder.to_dataframe()
