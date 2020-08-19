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
    COUNT_MISMATCHES_AFTER_VARIANT,
    PROTEIN_SEQUENCE_LENGTH,
    MAX_PROTEIN_SEQUENCES_PER_VARIANT,
    MIN_VARIANT_SEQUENCE_COVERAGE,
    VARIANT_SEQUENCE_ASSEMBLY,
    MIN_VARIANT_SEQUENCE_ASSEMBLY_OVERLAP_SIZE,
)

from .genetic_code import translate_cdna
from .protein_sequence_helpers import (
    sort_protein_sequences,
    group_equivalent_translations
)
from .reference_context_helpers import reference_contexts_for_variant
from .translation import Translation
from .translation_helpers import find_mutant_amino_acid_interval
from .value_object import ValueObject
from .variant_sequence_creator import VariantSequenceCreator
from .variant_orf_helpers import match_variant_sequence_to_reference_context

from .logging import get_logger

logger = get_logger(__name__)


class ProteinSequenceCreator(ValueObject):
    """
    Creates ProteinSequence objects for each variant by translating
    cDNA into one or more Translation objects and then grouping them
    by identical amino acid sequences. Each Translation comes from the
    combination of a variant cDNA sequence with the reading frames of
    annotated reference transcripts.
    """

    def __init__(
            self,
            protein_sequence_length=PROTEIN_SEQUENCE_LENGTH,
            min_variant_sequence_coverage=MIN_VARIANT_SEQUENCE_COVERAGE,
            min_transcript_prefix_length=MIN_TRANSCRIPT_PREFIX_LENGTH,
            max_transcript_mismatches=MAX_REFERENCE_TRANSCRIPT_MISMATCHES,
            count_mismatches_after_variant=COUNT_MISMATCHES_AFTER_VARIANT,
            max_protein_sequences_per_variant=MAX_PROTEIN_SEQUENCES_PER_VARIANT,
            variant_sequence_assembly=VARIANT_SEQUENCE_ASSEMBLY,
            min_assembly_overlap_size=MIN_VARIANT_SEQUENCE_ASSEMBLY_OVERLAP_SIZE):
        """
        protein_sequence_length : int
            Try to translate protein sequences of this length, though sometimes
            we'll have to return something shorter (depending on the RNAseq data,
            and presence of stop codons).

        min_variant_sequence_coverage : int
            Trim variant sequences to positions supported by at least this number
            of RNA reads.

        min_transcript_prefix_length : int
            Minimum number of bases we need to try matching between the reference
            context and variant sequence.

        max_transcript_mismatches : int
            Don't try to determine the reading frame for a transcript if more
            than this number of bases differ.

        count_mismatches_after_variant : bool
            Include mismatches after the variant locus in the count compared
            against max_transcript_mismatches.

        max_protein_sequences_per_variant : int
            Number of protein sequences to return for each ProteinSequence

        variant_sequence_assembly : bool
            If True, then assemble variant cDNA sequences based on overlap of
            RNA reads. If False, then variant cDNA sequences must be fully spanned
            and contained within RNA reads.

        min_assembly_overlap_size : int
            Minimum number of nucleotides that two reads need to overlap before they
            can be merged into a single coding sequence.
        """
        self.protein_sequence_length = protein_sequence_length
        self.min_variant_sequence_coverage = min_variant_sequence_coverage
        self.min_transcript_prefix_length = min_transcript_prefix_length
        self.max_transcript_mismatches = max_transcript_mismatches
        self.count_mismatches_after_variant = count_mismatches_after_variant
        self.variant_sequence_assembly = variant_sequence_assembly
        self.min_assembly_overlap_size = min_assembly_overlap_size

        # Adding an extra codon to the desired RNA sequence length in case we
        # need to clip nucleotides at the start/end of the sequence
        self._cdna_sequence_length = (self.protein_sequence_length + 1) * 3

        self._variant_sequence_creator = VariantSequenceCreator(
            min_variant_sequence_coverage=self.min_variant_sequence_coverage,
            preferred_sequence_length=self._cdna_sequence_length,
            variant_sequence_assembly=self.variant_sequence_assembly,
            min_assembly_overlap_size=self.min_assembly_overlap_size)

        self.max_protein_sequences_per_variant = max_protein_sequences_per_variant

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

        logger.info(
            "Full mutant cDNA sequence: %s (len=%d)",
            variant_sequence.sequence,
            len(variant_sequence))
        variant_orf = match_variant_sequence_to_reference_context(
            variant_sequence,
            reference_context,
            min_transcript_prefix_length=self.min_transcript_prefix_length,
            max_transcript_mismatches=self.max_transcript_mismatches,
            count_mismatches_after_variant=self.count_mismatches_after_variant)

        if variant_orf is None:
            logger.info("Unable to determine reading frame for %s", variant_sequence)
            return None

        cdna_sequence = variant_orf.cdna_sequence
        cdna_codon_offset = variant_orf.offset_to_first_complete_codon
        logger.info(
            "Untrimmed cDNA sequence: %s, offset to first codon = %d, len=%d",
            cdna_sequence,
            cdna_codon_offset,
            len(cdna_sequence))
        # get the offsets into the cDNA sequence which pick out the variant nucleotides
        cdna_variant_start_offset = variant_orf.variant_cdna_interval_start
        cdna_variant_end_offset = variant_orf.variant_cdna_interval_end

        in_frame_cdna_sequence = cdna_sequence[cdna_codon_offset:]
        logger.info("Translating '%s' (len=%d, expected AA length=%d)",
            in_frame_cdna_sequence,
            len(in_frame_cdna_sequence),
            len(in_frame_cdna_sequence) // 3)
        # TODO:
        #  determine if the first codon is the start codon of a
        #  transcript, for now any of the unusual start codons like CTG
        #  will translate to leucine instead of methionine
        amino_acids, ends_with_stop_codon = translate_cdna(
            in_frame_cdna_sequence,
            first_codon_is_start=False,
            mitochondrial=reference_context.mitochondrial)
        logger.info("Translated amino acids: %s, ends_with_stop=%s, len=%d" % (
            amino_acids,
            ends_with_stop_codon,
            len(amino_acids)))
        mutation_start_idx, mutation_end_idx, frameshift = \
            find_mutant_amino_acid_interval(
                cdna_sequence=cdna_sequence,
                cdna_first_codon_offset=cdna_codon_offset,
                cdna_variant_start_offset=cdna_variant_start_offset,
                cdna_variant_end_offset=cdna_variant_end_offset,
                n_ref=len(reference_context.sequence_at_variant_locus),
                n_amino_acids=len(amino_acids))

        if self.protein_sequence_length:
            if len(amino_acids) > self.protein_sequence_length:
                # if the protein is too long then shorten it, which implies
                # we're no longer stopping due to a stop codon and that the variant
                # amino acids might need a new stop index
                amino_acids = amino_acids[:self.protein_sequence_length]
                mutation_end_idx = min(
                    mutation_end_idx,
                    self.protein_sequence_length)
                ends_with_stop_codon = False

        if mutation_end_idx == mutation_start_idx:
            # a deletion only counts as mutated if the amino acids on
            # either side are in the sequence
            contains_mutation = (0 < mutation_start_idx < len(amino_acids))
        else:
            contains_mutation = len(amino_acids) > mutation_start_idx

        translation = Translation(
            amino_acids=amino_acids,
            contains_mutation=contains_mutation,
            frameshift=frameshift,
            ends_with_stop_codon=ends_with_stop_codon,
            mutation_start_idx=mutation_start_idx,
            mutation_end_idx=mutation_end_idx,
            untrimmed_variant_sequence=variant_sequence,
            reference_context=reference_context,
            variant_orf=variant_orf)

        logger.info(
            ("Translation from:"
             "\n-- cDNA = %s"
             "\n-- context = %s"
             "\n-- translation = %s") % (
                variant_sequence,
                reference_context,
                translation))

        return translation

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

        reference_contexts = reference_contexts_for_variant(
            variant,
            context_size=self._cdna_sequence_length,
            transcript_id_whitelist=transcript_id_whitelist)

        if len(reference_contexts) == 0:
            logger.info("Could not determine reference context for variant %s", variant)
            return []

        variant_sequences = self._variant_sequence_creator.reads_to_variant_sequences(
            variant=variant,
            reads=variant_reads)

        if not variant_sequences:
            logger.info("No spanning cDNA sequences for variant %s", variant)
            return []

        return self.all_pairs_translations(
            variant_sequences=variant_sequences,
            reference_contexts=reference_contexts)

    def translate_variants(
                self,
                variants_with_read_evidence_generator,
                transcript_id_whitelist=None):
            """
            Translates each coding variant in a collection to one or more protein
            fragment sequences (if the variant is not filtered and its spanning RNA
            sequences can be given a reading frame).

            Parameters
            ----------
            variants_with_read_evidence_generator : sequence or generator
                Each item of this sequence should be a pair containing a varcode.Variant
                and a ReadEvidence object

            transcript_id_whitelist : set, optional
                If given, expected to be a set of transcript IDs which we should use
                for determining the reading frame around a variant. If omitted, then
                try to use all overlapping reference transcripts.

            Yields pairs of a Variant and a sequence of all its candidate
            Translation objects.
            """
            for variant, read_evidence in variants_with_read_evidence_generator:
                translations = self.translate_variant_reads(
                    variant=variant,
                    variant_reads=read_evidence.alt_reads,
                    transcript_id_whitelist=transcript_id_whitelist)
                yield variant, translations

    def sorted_protein_sequences_for_variant(
            self,
            variant,
            read_evidence,
            transcript_id_whitelist=None):
        """"
        Translates a coding variant and its overlapping RNA reads into Translation
        objects, which are aggregated into ProteinSequence objects by their
        amino acid sequence (when they have equivalent coding sequences).

        Parameters
        ----------
        variant : varcode.Variant

        read_evidence : ReadEvidence object

        transcript_id_whitelist : set, optional
            If given, expected to be a set of transcript IDs which we should use
            for determining the reading frame around a variant. If omitted, then
            try to use all overlapping reference transcripts.

        Returns a list of ProteinSequence objects
        """
        translations = self.translate_variant_reads(
            variant=variant,
            variant_reads=read_evidence.alt_reads,
            transcript_id_whitelist=transcript_id_whitelist)

        # group distinct cDNA translations into ProteinSequence objects
        # by their amino acid sequence
        protein_sequences = group_equivalent_translations(translations)

        # sort protein sequences before returning the top results
        protein_sequences = sort_protein_sequences(protein_sequences)
        return protein_sequences

    def protein_sequences_from_read_evidence_generator(
            self,
            read_evidence_generator,
            transcript_id_whitelist=None):
        """

        Parameters
        ----------
        read_evidence_generator : generator of (varcode.Variant, ReadEvidence)
            Generator which yields sequence of Variant objects paired with
            their corresponding ReadEvidence

        transcript_id_whitelist : set of str or None
            Which transcripts should be considered when predicting DNA-only
            coding effects of mutations and also when trying to establish a
            reading frame for identified cDNA sequences.

        Generates sequence of (varcode.Variant, ProteinSequence list) pairs.
        """
        for variant, read_evidence in read_evidence_generator:
            protein_sequences = \
                self.sorted_protein_sequences_for_variant(
                    variant=variant,
                    read_evidence=read_evidence,
                    transcript_id_whitelist=transcript_id_whitelist)
            if self.max_protein_sequences_per_variant:
                protein_sequences = protein_sequences[:self.max_protein_sequences_per_variant]
            yield variant, protein_sequences
