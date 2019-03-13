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
from .genetic_code import translate_cdna
from .logging import get_logger
from .variant_sequence_creator import VariantSequenceCreator
from .reference_context import reference_contexts_for_variant
from .translation import Translation
from .translation_helpers import find_mutant_amino_acid_interval
from .variant_orf_helpers import match_variant_sequence_to_reference_context

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

        # Adding an extra codon to the desired RNA sequence length in case we
        # need to clip nucleotides at the start/end of the sequence
        self._cdna_sequence_length =  (self.protein_sequence_length + 1) * 3

        self._variant_sequence_creator = VariantSequenceCreator(
            min_variant_sequence_coverage=self.min_variant_sequence_coverage,
            preferred_sequence_length=self._cdna_sequence_length,
            variant_sequence_assembly=self.variant_sequence_assembly,
            min_assembly_overlap_size=self.min_assembly_overlap_size)

    def translation_from_variant_sequence_and_reference_context(
                self,
                variant_sequence,
                reference_context):
            """
            Attempt to translate a single VariantSequence using the reading frame
            from a single ReferenceContext.

            Parameters
            ----------
            variant_sequence : VariantSequence

            reference_context : ReferenceContext

            Returns either a Translation object or None if the number of
            mismatches between the RNA and reference transcript sequences exceeds
            given threshold.
            """
            variant_sequence_in_reading_frame = match_variant_sequence_to_reference_context(
                variant_sequence,
                reference_context,
                min_transcript_prefix_length=self.min_transcript_prefix_length,
                max_transcript_mismatches=self.max_transcript_mismatches,
                include_mismatches_after_variant=self.include_mismatches_after_variant)

            if variant_sequence_in_reading_frame is None:
                logger.info("Unable to determine reading frame for %s", variant_sequence)
                return None

            cdna_sequence = variant_sequence_in_reading_frame.cdna_sequence
            cdna_codon_offset = variant_sequence_in_reading_frame.offset_to_first_complete_codon

            # get the offsets into the cDNA sequence which pick out the variant nucleotides
            cdna_variant_start_offset = variant_sequence_in_reading_frame.variant_cdna_interval_start
            cdna_variant_end_offset = variant_sequence_in_reading_frame.variant_cdna_interval_end

            # TODO: determine if the first codon is the start codon of a
            # transcript, for now any of the unusual start codons like CTG
            # will translate to leucine instead of methionine.
            variant_amino_acids, ends_with_stop_codon = translate_cdna(
                cdna_sequence[cdna_codon_offset:],
                first_codon_is_start=False,
                mitochondrial=reference_context.mitochondrial)

            variant_aa_interval_start, variant_aa_interval_end, frameshift = \
                find_mutant_amino_acid_interval(
                    cdna_sequence=cdna_sequence,
                    cdna_first_codon_offset=cdna_codon_offset,
                    cdna_variant_start_offset=cdna_variant_start_offset,
                    cdna_variant_end_offset=cdna_variant_end_offset,
                    n_ref=len(reference_context.sequence_at_variant_locus),
                    n_amino_acids=len(variant_amino_acids))

            if self.protein_sequence_length:
                if len(variant_amino_acids) > self.protein_sequence_length:
                    if self.protein_sequence_length <= variant_aa_interval_start:
                        logger.warn(
                            ("Truncating amino acid sequence %s "
                             "to only %d elements loses all variant residues"),
                            variant_amino_acids,
                            self.protein_sequence_length)
                        return None
                    else:
                        # if the protein is too long then shorten it, which implies
                        # we're no longer stopping due to a stop codon and that the variant
                        # amino acids might need a new stop index
                        variant_amino_acids = variant_amino_acids[:self.protein_sequence_length]
                        variant_aa_interval_end = min(
                            variant_aa_interval_end,
                            self.protein_sequence_length)
                        ends_with_stop_codon = False

            return Translation(
                amino_acids=variant_amino_acids,
                frameshift=frameshift,
                ends_with_stop_codon=ends_with_stop_codon,
                variant_aa_interval_start=variant_aa_interval_start,
                variant_aa_interval_end=variant_aa_interval_end,
                untrimmed_variant_sequence=variant_sequence,
                reference_context=reference_context,
                variant_sequence_in_reading_frame=variant_sequence_in_reading_frame)

    def all_pairs_translations(
            self,
            variant_sequences,
            reference_contexts):
        """
        Given all a list of VariantSequence objects for a particular variant
        and all the ReferenceContext objects for that locus, attempt to
        translate all pairs of sequences and reference contexts.

        Parameters
        ----------
        variant_sequences : list of VariantSequence

        reference_contexts : list of ReferenceContext

        Return list of Translation objects.
        """
        translations = []
        for reference_context in reference_contexts:
            for variant_sequence in variant_sequences:
                translation = self.translation_from_variant_sequence_and_reference_context(
                    variant_sequence=variant_sequence,
                    reference_context=reference_context)
                if translation is not None:
                    translations.append(translation)
        return translations

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

        Returns list of Translation objects
        """
        if len(variant_reads) == 0:
            logger.info("No supporting reads for variant %s", variant)
            return []

        variant_sequences = self._variant_sequence_creator.reads_to_variant_sequences(
            variant=variant,
            reads=variant_reads)

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

        return self.all_pairs_translations(
            variant_sequences=variant_sequences,
            reference_contexts=reference_contexts)

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
                translations = self.translate_variant_reads(
                    variant=variant,
                    variant_reads=variant_reads,
                    transcript_id_whitelist=transcript_id_whitelist)
                yield variant, translations
