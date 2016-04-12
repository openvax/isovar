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
        # number of unique read names from all the VariantSequence objects
        # from each translation
        "number_supporting_reads",
        # total number of reads at the locus which contained the variant
        # nucleotides, even if they supported other phased sequences
        "total_variant_reads",
        # how many reference transcripts were used to establish the
        # reading frame for this protein sequence
        "number_supporting_transcripts",
        # how many reference transcripts overlap the variant locus?
        "total_overlapping_transcripts"
        # name of gene of the reference transcripts used in Translation
        # objects
        "gene",
    ))

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
        # dictionary mapping TranslationKey object to a list of Translations
        # with the same protein sequence
        equivalent_translations_dict = defaultdict(list)
        for translation in translations:
            translation_key_values = {
                name: getattr(translation, name)
                for name in TranslationKey._fields
            }
            key = TranslationKey(**translation_key_values)
            equivalent_translations_dict[key].append(translation)

        protein_sequences = []
        for (key, equivalent_translations) in equivalent_translations_dict.items():
            # first fill in all the ProteinSequence fields which are shared
            # with TranslationKey
            protein_sequence_values = {
                name: getattr(key, name)
                for name in TranslationKey._fields
            }
            protein_sequence_values["translations"] = equivalent_translations
            # TODO: finish this
            protein_sequence_values["number_supporting_reads"] = None
            protein_sequence_values["total_variant_reads"] = None
            protein_sequence_values["number_supporting_transcripts"] = None
            protein_sequence_values["total_overlapping_transcripts"] = None
            protein_sequence_values["gene"] = None
            protein_sequences.append(ProteinSequence(**protein_sequence_values))
            # TODO: sort protein sequences
        yield variant, protein_sequences[:max_protein_sequences_per_variant]

def variants_to_protein_sequences_dataframe(*args, **kwargs):
    """
    Given a collection of variants and a SAM/BAM file of overlapping reads,
    returns a DataFrame with a row for each protein sequence.

    Takes the same parameters as variants_to_protein_sequences.
    """
    df_builder = DataFrameBuilder(
        key_by_variant=True,
        column_names=ProteinSequence._fields,
        exclude=["translations"])
    for (variant, protein_sequences) in variants_to_protein_sequences(*args, **kwargs):
        for protein_sequence in protein_sequences:
            df_builder.add(variant, protein_sequence)
    return df_builder.to_dataframe()
