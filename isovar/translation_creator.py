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

from .default_parameters import (
    MIN_TRANSCRIPT_PREFIX_LENGTH,
    MAX_REFERENCE_TRANSCRIPT_MISMATCHES,
    INCLUDE_MISMATCHES_AFTER_VARIANT,
    PROTEIN_SEQUENCE_LENGTH,
    MIN_ALT_RNA_READS,
    MIN_VARIANT_SEQUENCE_COVERAGE,
    VARIANT_SEQUENCE_ASSEMBLY,
    MIN_VARIANT_SEQUENCE_ASSEMBLY_OVERLAP_SIZE
)
from .logging import get_logger

logger = get_logger(__name__)

class TranslationCreator(object):
    """
    TranslationCreator is used to combine variant cDNA sequences from a BAM
    file with the reading frames of annotated reference transcripts to create
    candidate translations.
    """
    def __init__(
            self,
            protein_sequence_length=PROTEIN_SEQUENCE_LENGTH,
            min_alt_rna_reads=MIN_ALT_RNA_READS,
            min_variant_sequence_coverage=MIN_VARIANT_SEQUENCE_COVERAGE,
            min_transcript_prefix_length=MIN_TRANSCRIPT_PREFIX_LENGTH,
            max_transcript_mismatches=MAX_REFERENCE_TRANSCRIPT_MISMATCHES,
            include_mismatches_after_variant=INCLUDE_MISMATCHES_AFTER_VARIANT,
            variant_sequence_assembly=VARIANT_SEQUENCE_ASSEMBLY,
            min_assembly_overlap_size=MIN_VARIANT_SEQUENCE_ASSEMBLY_OVERLAP_SIZE):
        """
        Parameters
        ----------
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
        """
        self.protein_sequence_length = protein_sequence_length
        self.min_alt_rna_reads = min_alt_rna_reads
        self.min_variant_sequence_coverage = min_variant_sequence_coverage
        self.min_transcript_prefix_length = min_transcript_prefix_length
        self.max_transcript_mismatches = max_transcript_mismatches
        self.include_mismatches_after_variant = include_mismatches_after_variant
        self.variant_sequence_assembly = variant_sequence_assembly
        self.min_assembly_overlap_size = min_assembly_overlap_size

    def translate_variant_reads(
            self,
            variant,
            variant_reads,
            transcript_id_whitelist=None):
        """
        Given a variant and its associated alt reads, construct variant sequences
        and translate them into Translation objects.

        Returns 0 or more Translation objects.

        Parameters
        ----------
        variant : varcode.Variant

        variant_reads : sequence or generator
            AlleleRead objects supporting the variant

        transcript_id_whitelist : set, optional
            If given, expected to be a set of transcript IDs which we should use
            for determining the reading frame around a variant. If omitted, then
            try to use all overlapping reference transcripts.
        """
        if len(variant_reads) == 0:
            logger.info("No supporting reads for variant %s", variant)
            return []

        # Adding an extra codon to the desired RNA sequence length in case we
        # need to clip nucleotides at the start/end of the sequence
        cdna_sequence_length = (self.protein_sequence_length + 1) * 3

        variant_sequences = reads_to_variant_sequences(
            variant=variant,
            reads=variant_reads,
            preferred_sequence_length=cdna_sequence_length,
            min_alt_rna_reads=self.min_alt_rna_reads,
            min_variant_sequence_coverage=self.min_variant_sequence_coverage,
            variant_sequence_assembly=self.variant_sequence_assembly)

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
            reference_contexts=reference_contexts))
    def translate_variants(
                self,
                variants_with_supporting_reads,
                transcript_id_whitelist=None):
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

            Yields pairs of a Variant and a sequence of all its candidate
            Translation objects.
            """
            for variant, variant_reads in variants_with_supporting_reads:
                translations = translate_variant_reads(
                    variant=variant,
                    variant_reads=variant_reads,
                    transcript_id_whitelist=transcript_id_whitelist,
                    min_alt_rna_reads=min_alt_rna_reads,
                    min_variant_sequence_coverage=min_variant_sequence_coverage,
                    min_transcript_prefix_length=min_transcript_prefix_length,
                    max_transcript_mismatches=max_transcript_mismatches,
                    include_mismatches_after_variant=include_mismatches_after_variant,
                    variant_sequence_assembly=variant_sequence_assembly)
                yield variant, translations
