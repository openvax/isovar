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

"""
This module combines variant cDNA sequences collected from a BAM file with
the reading frames of annotated reference transcripts to create candidate
translations.
"""


from __future__ import print_function, division, absolute_import

from .default_parameters import (
    MIN_TRANSCRIPT_PREFIX_LENGTH,
    MAX_REFERENCE_TRANSCRIPT_MISMATCHES,
    INCLUDE_MISMATCHES_AFTER_VARIANT,
    PROTEIN_SEQUENCE_LENGTH,
    MIN_ALT_RNA_READS,
    MIN_VARIANT_SEQUENCE_COVERAGE,
    VARIANT_SEQUENCE_ASSEMBLY,
)
from .reference_context import reference_contexts_for_variant
from .variant_sequence_helpers import reads_to_variant_sequences


def translation_generator(
        variant_sequences,
        reference_contexts,
        min_transcript_prefix_length,
        max_transcript_mismatches,
        include_mismatches_after_variant,
        protein_sequence_length=None):
    """
    Given all detected VariantSequence objects for a particular variant
    and all the ReferenceContext objects for that locus, translate
    multiple protein sequences, up to the number specified by the argument
    max_protein_sequences_per_variant.

    Parameters
    ----------
    variant_sequences : list of VariantSequence objects
        Variant sequences overlapping a single original variant

    reference_contexts : list of ReferenceContext objects
        Reference sequence contexts from the same variant as the variant_sequences

    min_transcript_prefix_length : int
        Minimum number of nucleotides before the variant to test whether
        our variant sequence can use the reading frame from a reference
        transcript.

    max_transcript_mismatches : int
        Maximum number of mismatches between coding sequence before variant
        and reference transcript we're considering for determing the reading
        frame.

    include_mismatches_after_variant : bool
        If true, mismatches occurring after the variant locus will also count
        toward max_transcript_mismatches filtering.

    protein_sequence_length : int, optional
        Truncate protein to be at most this long.

    Yields a sequence of Translation objects.
    """
    for reference_context in reference_contexts:
        for variant_sequence in variant_sequences:
            translation = Translation.from_variant_sequence_and_reference_context(
                variant_sequence=variant_sequence,
                reference_context=reference_context,
                min_transcript_prefix_length=min_transcript_prefix_length,
                max_transcript_mismatches=max_transcript_mismatches,
                include_mismatches_after_variant=include_mismatches_after_variant,
                protein_sequence_length=protein_sequence_length)
            if translation is not None:
                yield translation


def translate_variant_reads(
        variant,
        variant_reads,
        protein_sequence_length,
        transcript_id_whitelist=None,
        min_alt_rna_reads=MIN_ALT_RNA_READS,
        min_variant_sequence_coverage=MIN_VARIANT_SEQUENCE_COVERAGE,
        min_transcript_prefix_length=MIN_TRANSCRIPT_PREFIX_LENGTH,
        max_transcript_mismatches=MAX_REFERENCE_TRANSCRIPT_MISMATCHES,
        include_mismatches_after_variant=INCLUDE_MISMATCHES_AFTER_VARIANT,
        variant_sequence_assembly=VARIANT_SEQUENCE_ASSEMBLY):
    """
    Given a variant and its associated alt reads, construct variant sequences
    and translate them into Translation objects.

    Returns 0 or more Translation objects.

    Parameters
    ----------
    variant : varcode.Variant

    variant_reads : sequence or generator
        AlleleRead objects supporting the variant

    protein_sequence_length : int
        Try to translate protein sequences of this length, though sometimes
        we'll have to return something shorter (depending on the RNAseq data,
        and presence of stop codons).

    transcript_id_whitelist : set, optional
        If given, expected to be a set of transcript IDs which we should use
        for determining the reading frame around a variant. If omitted, then
        try to use all overlapping reference transcripts.

    min_alt_rna_reads : int
        Drop variant sequences from loci with fewer than this number of
        RNA reads supporting the alt allele.

    min_variant_sequence_coverage : int
        Trim variant sequences to nucleotides covered by at least this many
        reads.

    min_transcript_prefix_length : int
        Minimum number of bases we need to try matching between the reference
        context and variant sequence.

    max_transcript_mismatches : int
        Don't try to determine the reading frame for a transcript if more
        than this number of bases differ.

    include_mismatches_after_variant : bool
        Include mismatches after the variant locus in the count compared
        against max_transcript_mismatches.

    variant_sequence_assembly : bool
        Use overlap assembly to construct longer variant cDNA sequences.
    """
    if len(variant_reads) == 0:
        logger.info("No supporting reads for variant %s", variant)
        return []

    # Adding an extra codon to the desired RNA sequence length in case we
    # need to clip nucleotides at the start/end of the sequence
    cdna_sequence_length = (protein_sequence_length + 1) * 3

    variant_sequences = reads_to_variant_sequences(
        variant=variant,
        reads=variant_reads,
        preferred_sequence_length=cdna_sequence_length,
        min_alt_rna_reads=min_alt_rna_reads,
        min_variant_sequence_coverage=min_variant_sequence_coverage,
        variant_sequence_assembly=variant_sequence_assembly)

    if not variant_sequences:
        logger.info("No spanning cDNA sequences for variant %s", variant)
        return []

    # try translating the variant sequences from the same set of
    # ReferenceContext objects, which requires using the longest
    # context_size to be compatible with all of the sequences. Some
    # sequences maybe have fewer nucleotides than this before the variant
    # and will thus have to be trimmed.
    context_size = max(
        len(variant_sequence.prefix)
        for variant_sequence in variant_sequences)

    reference_contexts = reference_contexts_for_variant(
        variant,
        context_size=context_size,
        transcript_id_whitelist=transcript_id_whitelist)

    return list(translation_generator(
        variant_sequences=variant_sequences,
        reference_contexts=reference_contexts,
        min_transcript_prefix_length=min_transcript_prefix_length,
        max_transcript_mismatches=max_transcript_mismatches,
        include_mismatches_after_variant=include_mismatches_after_variant,
        protein_sequence_length=protein_sequence_length))


def translate_variants(
        variants_with_supporting_reads,
        transcript_id_whitelist=None,
        protein_sequence_length=PROTEIN_SEQUENCE_LENGTH,
        min_alt_rna_reads=MIN_ALT_RNA_READS,
        min_variant_sequence_coverage=MIN_VARIANT_SEQUENCE_COVERAGE,
        min_transcript_prefix_length=MIN_TRANSCRIPT_PREFIX_LENGTH,
        max_transcript_mismatches=MAX_REFERENCE_TRANSCRIPT_MISMATCHES,
        include_mismatches_after_variant=INCLUDE_MISMATCHES_AFTER_VARIANT,
        variant_sequence_assembly=VARIANT_SEQUENCE_ASSEMBLY):
    """
    Translates each coding variant in a collection to one or more protein
    fragment sequences (if the variant is not filtered and its spanning RNA
    sequences can be given a reading frame).

    Parameters
    ----------
    variants_with_supporting_reads : sequence or generator
        Each item of this sequence should be a pair containing a varcode.Variant
        and a list of AlleleRead objects supporting that variant.

    transcript_id_whitelist : set, optional
        If given, expected to be a set of transcript IDs which we should use
        for determining the reading frame around a variant. If omitted, then
        try to use all overlapping reference transcripts.

    protein_sequence_length : int
        Try to translate protein sequences of this length, though sometimes
        we'll have to return something shorter (depending on the RNAseq data,
        and presence of stop codons).

    min_alt_rna_reads : int
        Drop variant sequences from loci with fewer than this number of
        RNA reads supporting the alt allele.

    min_variant_sequence_coverage : int
        Trim variant sequences to nucleotides covered by at least this many
        reads.

    min_transcript_prefix_length : int
        Minimum number of bases we need to try matching between the reference
        context and variant sequence.

    max_transcript_mismatches : int
        Don't try to determine the reading frame for a transcript if more
        than this number of bases differ.

    include_mismatches_after_variant : bool
        Include mismatches after the variant locus in the count compared
        against max_transcript_mismatches.

    variant_sequence_assembly : bool
        Use overlap assembly to construct longer variant cDNA sequences.

    Yields pairs of a Variant and a sequence of all its candidate
    Translation objects.
    """
    for variant, variant_reads in variants_with_supporting_reads:
        translations = translate_variant_reads(
            variant=variant,
            variant_reads=variant_reads,
            protein_sequence_length=protein_sequence_length,
            transcript_id_whitelist=transcript_id_whitelist,
            min_alt_rna_reads=min_alt_rna_reads,
            min_variant_sequence_coverage=min_variant_sequence_coverage,
            min_transcript_prefix_length=min_transcript_prefix_length,
            max_transcript_mismatches=max_transcript_mismatches,
            include_mismatches_after_variant=include_mismatches_after_variant,
            variant_sequence_assembly=variant_sequence_assembly)
        yield variant, translations
