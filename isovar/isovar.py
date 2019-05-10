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
    MIN_ALT_RNA_READS,
    MIN_ALT_RNA_FRAGMENTS,
    MIN_RNA_VAF,
    MIN_RATIO_ALT_TO_OTHER_NONREF_RNA_FRAGMENTS,
    MIN_VARIANT_SEQUENCE_COVERAGE,
    VARIANT_SEQUENCE_ASSEMBLY,
    MIN_VARIANT_SEQUENCE_ASSEMBLY_OVERLAP_SIZE,
    USE_SECONDARY_ALIGNMENTS,
    USE_DUPLICATE_READS,
    MIN_READ_MAPPING_QUALITY,
    USE_SOFT_CLIPPED_BASES
)

from .protein_sequence_helpers import sort_protein_sequences, collapse_translations
from .protein_sequence_creator import ProteinSequenceCreator
from .read_collector import ReadCollector
from .logging import get_logger
from .isovar_result import IsovarResult

logger = get_logger(__name__)

class Isovar(object):
    """
    This is the main entrypoint into the Isovar library, which collects
    RNA reads supporting variants and translates their coding sequence
    into amino acid sequences.
    """
    def __init__(
            self,
            read_collector=None,
            protein_sequence_creator=None,
            min_alt_rna_fragments=MIN_ALT_RNA_FRAGMENTS,
            min_alt_rna_reads=MIN_ALT_RNA_READS,
            min_rna_vaf=MIN_RNA_VAF,
            min_ratio_alt_to_other_rna_fragments=MIN_RATIO_ALT_TO_OTHER_NONREF_RNA_FRAGMENTS,
            max_other_vaf=MAX_OTHER_VAF):
        """
        protein_sequence_length : int
            Try to translate protein sequences of this length, though sometimes
            we'll have to return something shorter (depending on the RNAseq data,
            and presence of stop codons).

        min_alt_rna_reads : int
            Drop variant sequences at loci with fewer than this number of reads
            supporting the alt allele.

        min_ratio_alt_to_other_nonref_reads : float
            Drop variant sequences at loci where there is support for third/fourth
            alleles in the RNA and the count for the alt allele is not at least
            this number greater than the sum of other non-ref alleles.

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

        use_secondary_alignments : bool
            Use a read even when it's not the primary alignment at a locus

        use_duplicate_reads : bool
            Use a read even if it's been marked as a duplicate

        min_mapping_quality : int
            Minimum MAPQ (mapping quality) to use a read

        use_soft_clipped_bases : bool
            Include soft-clipped positions on a read which were ignored by the aligner
        """
        if read_collector is None:
            read_collector = ReadCollector()

        self.read_collector = read_collector

        if protein_sequence_creator is None:
            self.protein_sequence_creator = ProteinSequenceCreator()

        self.min_alt_rna_fragments = min_alt_rna_fragments
        self.min_rna_vaf = min_rna_vaf
        self.min_alt_rna_reads = min_alt_rna_reads
        self.min_ratio_alt_to_other_rna_fragments = min_ratio_alt_to_other_rna_fragments
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
    
    def process_variants(
            self,
            variants,
            aligned_rna_reads,
            transcript_id_whitelist=None):
        """

        Parameters
        ----------
        variants : varcode.VariantCollection
            Somatic variants

        aligned_rna_reads : pysam.AlignmentFile
            Aligned tumor RNA reads

        transcript_id_whitelist : set of str or None
            Which transcripts should be considered when predicting DNA-only
            coding effects of mutations and also when trying to establish a
            reading frame for identified cDNA sequences.

        Generator of IsovarResult objects, one for each variant.
        The `protein_sequences` field of the IsovarVar result will be empty
        if no sequences could be determined.
        """

        # create generator which returns (Variant, GroupedAlleleReads) pairs
        variant_and_read_gen = \
            self._read_collector.grouped_allele_reads_overlapping_variants(
               variants=variants,
               alignment_file=aligned_rna_reads)

        for variant, grouped_allele_reads in  variant_and_read_gen:
            # generate protein sequences by assembling variant reads
            protein_sequences = \
                self.sorted_protein_sequences_for_variant(
                    variant=variant,
                    grouped_allele_reads=grouped_allele_reads,
                    transcript_id_whitelist=transcript_id_whitelist)

            # if we're only keeping the top-k, get rid of the rest
            if self.max_protein_sequences_per_variant:
                protein_sequences = protein_sequences[:self.max_protein_sequences_per_variant]

            yield IsovarResult(
                variant=variant,
                transcript_id_whitelist=transcript_id_whitelist,
                grouped_allele_reads=grouped_allele_reads,
                protein_sequences=protein_sequences)