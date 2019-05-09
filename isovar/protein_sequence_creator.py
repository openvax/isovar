# Copyright (c) 2019. Mount Sinai School of Medicine
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
    MAX_PROTEIN_SEQUENCES_PER_VARIANT,

    MIN_VARIANT_SEQUENCE_COVERAGE,
    VARIANT_SEQUENCE_ASSEMBLY,
    MIN_VARIANT_SEQUENCE_ASSEMBLY_OVERLAP_SIZE,
    USE_SECONDARY_ALIGNMENTS,
    USE_DUPLICATE_READS,
    MIN_READ_MAPPING_QUALITY,
    USE_SOFT_CLIPPED_BASES
)

from .protein_sequence_helpers import sort_protein_sequences, collapse_translations
from .translation_creator import TranslationCreator
from .logging import get_logger

from .read_collector import ReadCollector

logger = get_logger(__name__)


class ProteinSequenceCreator(TranslationCreator):
    """
    Creates ProteinSequence objects for each variant by translating
    cDNA into one or more Translation objects and then grouping them
    by identical amino acid sequences.
    """

    def __init__(
            self,
            protein_sequence_length=PROTEIN_SEQUENCE_LENGTH,
            min_variant_sequence_coverage=MIN_VARIANT_SEQUENCE_COVERAGE,
            min_transcript_prefix_length=MIN_TRANSCRIPT_PREFIX_LENGTH,
            max_transcript_mismatches=MAX_REFERENCE_TRANSCRIPT_MISMATCHES,
            include_mismatches_after_variant=INCLUDE_MISMATCHES_AFTER_VARIANT,
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

        include_mismatches_after_variant : bool
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
        TranslationCreator.__init__(
            self,
            protein_sequence_length=protein_sequence_length,
            min_variant_sequence_coverage=min_variant_sequence_coverage,
            min_transcript_prefix_length=min_transcript_prefix_length,
            max_transcript_mismatches=max_transcript_mismatches,
            include_mismatches_after_variant=include_mismatches_after_variant,
            variant_sequence_assembly=variant_sequence_assembly,
            min_assembly_overlap_size=min_assembly_overlap_size)

        self.max_protein_sequences_per_variant = max_protein_sequences_per_variant

    def sorted_protein_sequences_for_variant(
            self,
            variant,
            grouped_allele_reads,
            transcript_id_whitelist=None):
        """"
        Translates a coding variant and its overlapping RNA reads into Translation
        objects, which are aggregated into ProteinSequence objects by their
        amino acid sequence (when they have equivalent coding sequences).

        Parameters
        ----------
        variant : varcode.Variant

        grouped_allele_reads : GroupedAlleleReads object

        transcript_id_whitelist : set, optional
            If given, expected to be a set of transcript IDs which we should use
            for determining the reading frame around a variant. If omitted, then
            try to use all overlapping reference transcripts.

        Returns a list of ProteinSequence objects
        """
        translations = self.translate_variant_reads(
            variant=variant,
            variant_reads=grouped_allele_reads.alt_reads,
            transcript_id_whitelist=transcript_id_whitelist)

        # group distinct cDNA translations into ProteinSequence objects
        # by their amino acid sequence
        protein_sequences = collapse_translations(translations)

        # sort protein sequences before returning the top results
        protein_sequences = sort_protein_sequences(protein_sequences)
        return protein_sequences

    def variant_and_protein_sequences_generator(
            self,
            variant_and_reads_generator,
            transcript_id_whitelist=None):
        """

        Parameters
        ----------
        variant_and_reads_generator : generator of (varcode.Variant, GroupedAlleleReads)
            Generator which yields sequence of Variant objects paired with
            their corresponding GroupedAlleleReads

        transcript_id_whitelist : set of str or None
            Which transcripts should be considered when predicting DNA-only
            coding effects of mutations and also when trying to establish a
            reading frame for identified cDNA sequences.

        Generates sequence of (varcode.Variant, ProteinSequence list) pairs.
        """
        for variant, grouped_allele_reads in variant_and_reads_generator:
            protein_sequences = \
                self.sorted_protein_sequences_for_variant(
                    variant=variant,
                    grouped_allele_reads=grouped_allele_reads,
                    transcript_id_whitelist=transcript_id_whitelist)
            if self.max_protein_sequences_per_variant:
                protein_sequences = protein_sequences[:self.max_protein_sequences_per_variant]
            yield variant, protein_sequences
