# Copyright (c) 2016-2019. Mount Sinai School of Medicine
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

from .allele_counts import AlleleCount, count_alleles_at_variant_locus
from .allele_read import AlleleRead
from .common import list_to_string
from .dataframe_builder import DataFrameBuilder
from .locus_read import LocusRead, locus_read_generator
from .protein_sequence import ProteinSequence
from .reference_context import reference_contexts_for_variants, ReferenceContext
from .translation import Translation
from .variant_sequence import VariantSequence


def dataframe_from_generator(
        element_class,
        variant_and_elements_generator,
        **kwargs):
    builder = DataFrameBuilder(element_class, **kwargs)
    for variant, elements in variant_and_elements_generator:
        builder.add_many(variant, elements)
    return builder.to_dataframe()


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


def allele_counts_dataframe(variant_and_allele_reads_generator):
    """
    Creates a DataFrame containing number of reads supporting the
    ref vs. alt alleles for each variant.
    """
    df_builder = DataFrameBuilder(
        AlleleCount,
        extra_column_fns={
            "gene": lambda variant, _: ";".join(variant.gene_names),
        })
    for variant, allele_reads in variant_and_allele_reads_generator:
        counts = count_alleles_at_variant_locus(variant, allele_reads)
        df_builder.add(variant, counts)
    return df_builder.to_dataframe()


def allele_reads_to_dataframe(variants_and_allele_reads):
    """
    Parameters
    ----------
    variants_and_allele_reads : sequence
        List or generator of pairs whose first element is a Variant and
        whose second element is a sequence of AlleleRead objects.
    """
    df_builder = DataFrameBuilder(
        AlleleRead,
        extra_column_fns={
            "gene": lambda variant, _: ";".join(variant.gene_names),
        })
    for variant, allele_reads in variants_and_allele_reads:
        df_builder.add_many(variant, allele_reads)
    return df_builder.to_dataframe()


def locus_reads_dataframe(*args, **kwargs):
    """
    Traverse a BAM file to find all the reads overlapping a specified locus.

    Parameters are the same as those for read_locus_generator.
    """
    df_builder = DataFrameBuilder(
        LocusRead,
        variant_columns=False,
        converters={
            "reference_positions": list_to_string,
            "quality_scores": list_to_string,
        })
    for locus_read in locus_read_generator(*args, **kwargs):
        df_builder.add(variant=None, element=locus_read)
    return df_builder.to_dataframe()


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


def variant_sequences_generator_to_dataframe(variant_sequences_generator):
    """
    Creates a dataframe from a generator which yields
    (Variant, [VariantSequence]) pairs.

    Returns pandas.DataFrame
    """
    # TODO: Change VariantSequence.alt to VariantSequence.alt_nucleotides
    # or something else that doesn't clash with a variant's `alt` field
    return dataframe_from_generator(
        VariantSequence,
        variant_sequences_generator,
        rename_dict={"alt": "allele"},
        extra_column_fns={
            "gene": lambda variant, _: ";".join(variant.gene_names),
        })


def translations_generator_to_dataframe(translations_generator):
    """
    Given a generator of (Variant, [Translation]) pairs,
    returns a DataFrame of translated protein fragments with columns
    for each field of a Translation object (and chr/pos/ref/alt per variant).
    """
    return dataframe_from_generator(
        element_class=Translation,
        variant_and_elements_generator=translations_generator,
        exclude=[],
        converters={
            "untrimmed_variant_sequence": lambda vs: vs.sequence,
            "variant_sequence_in_reading_frame": (
                lambda vs: vs.in_frame_cdna_sequence),
            "reference_context": (
                lambda rc: ";".join([
                    transcript.name for
                    transcript in rc.transcripts]))
        },
        extra_column_fns={
            "untrimmed_variant_sequence_read_count": (
                lambda _, t: len(t.untrimmed_variant_sequence.reads)),
        })
