# Copyright (c) 2018. Mount Sinai School of Medicine
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
    VARIANT_SEQUENCE_ASSEMBLY
)
from .protein_sequence import ProteinSequence
from .protein_sequence_helpers import sort_protein_sequences
from .common import groupby
from .translation import translate_variant_reads, Translation
from .logging import get_logger
from .variant_support import gather_variant_support

logger = get_logger(__name__)


class ProteinSequenceCreator(object):
    def __init__(
            self,
            transcript_id_whitelist=None,
            protein_sequence_length=PROTEIN_SEQUENCE_LENGTH,
            min_alt_rna_fragments=MIN_ALT_RNA_FRAGMENTS,
            min_rna_vaf=MIN_RNA_VAF,
            min_alt_rna_reads=MIN_ALT_RNA_READS,
            min_ratio_alt_to_other_rna_fragments=MIN_RATIO_ALT_TO_OTHER_NONREF_RNA_FRAGMENTS,
            min_variant_sequence_coverage=MIN_VARIANT_SEQUENCE_COVERAGE,
            min_transcript_prefix_length=MIN_TRANSCRIPT_PREFIX_LENGTH,
            max_transcript_mismatches=MAX_REFERENCE_TRANSCRIPT_MISMATCHES,
            include_mismatches_after_variant=INCLUDE_MISMATCHES_AFTER_VARIANT,
            max_protein_sequences_per_variant=MAX_PROTEIN_SEQUENCES_PER_VARIANT,
            variant_sequence_assembly=VARIANT_SEQUENCE_ASSEMBLY):
        """
        transcript_id_whitelist : set, optional
            If given, expected to be a set of transcript IDs which we should use
            for determining the reading frame around a variant. If omitted, then
            try to use all overlapping reference transcripts.

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
        """
        self.transcript_id_whitelist = transcript_id_whitelist
        self.protein_sequence_length = protein_sequence_length
        self.min_alt_rna_fragments = min_alt_rna_fragments
        self.min_rna_vaf = min_rna_vaf
        self.min_alt_rna_reads = self.min_alt_rna_reads
        self.min_ratio_alt_to_other_rna_fragments = min_ratio_alt_to_other_rna_fragments
        self.min_variant_sequence_coverage = min_variant_sequence_coverage
        self.min_transcript_prefix_length = min_transcript_prefix_length
        self.max_transcript_mismatches = max_transcript_mismatches
        self.include_mismatches_after_variant = include_mismatches_after_variant
        self.max_protein_sequences_per_variant = max_protein_sequences_per_variant
        self.variant_sequence_assembly = variant_sequence_assembly

    def protein_sequences_from_variant_and_overlapping_reads(
            self,
            variant,
            overlapping_reads):
        """"
        Translates a coding variant and its overlapping RNA reads into Translation
        objects, which are aggregated into ProteinSequence objects by their
        amino acid sequence (when they have equivalent coding sequences).

        Parameters
        ----------
        variant : varcode.Variant

        overlapping_reads : list of AlleleReads

        Returns a list of ProteinSequence objects
        """
        variant_support = gather_variant_support(variant, overlapping_reads)

        translations = translate_variant_reads(
            variant=variant,
            variant_reads=variant_support.alt_reads,
            transcript_id_whitelist=self.transcript_id_whitelist,
            protein_sequence_length=self.protein_sequence_length,
            min_variant_sequence_coverage=self.min_variant_sequence_coverage,
            min_transcript_prefix_length=self.min_transcript_prefix_length,
            max_transcript_mismatches=self.max_transcript_mismatches,
            include_mismatches_after_variant=self.include_mismatches_after_variant,
            variant_sequence_assembly=self.variant_sequence_assembly)

        overlapping_transcript_ids = {
            t.id
            for t in variant.transcripts
            if t.is_protein_coding
        }
        protein_sequences = []
        for (key, equivalent_translations) in groupby(
                translations, key_fn=Translation.as_translation_key).items():

            # get the variant read names, transcript IDs and gene names for
            # protein sequence we're about to construct
            alt_reads_supporting_protein_sequence, group_transcript_ids, group_gene_names = \
                ProteinSequence._summarize_translations(equivalent_translations)

            logger.info(
                "%s: %s alt reads supporting protein sequence (gene names = %s)",
                key,
                len(alt_reads_supporting_protein_sequence),
                group_gene_names)

            protein_sequence = ProteinSequence.from_translation_key(
                translation_key=key,
                translations=equivalent_translations,
                overlapping_reads=overlapping_reads,
                alt_reads=variant_support.alt_reads,
                ref_reads=variant_support.ref_reads,
                alt_reads_supporting_protein_sequence=alt_reads_supporting_protein_sequence,
                transcripts_supporting_protein_sequence=group_transcript_ids,
                transcripts_overlapping_variant=overlapping_transcript_ids,
                gene=list(group_gene_names))
            logger.info("%s: protein sequence = %s" % (key, protein_sequence.amino_acids))
            protein_sequences.append(protein_sequence)
        # sort protein sequences before returning the top results
        protein_sequences = sort_protein_sequences(protein_sequences)
        return protein_sequences

    def reads_generator_to_protein_sequences_generator(
            self,
            variant_and_overlapping_reads_generator,
            transcript_id_whitelist=None):
        """"
        Translates each coding variant in a collection to one or more
        Translation objects, which are then aggregated into equivalent
        ProteinSequence objects.

        Parameters
        ----------
        variant_and_overlapping_reads_generator : generator
            Yields sequence of varcode.Variant objects paired with sequences
            of AlleleRead objects that support that variant.

        transcript_id_whitelist : set, optional
            If given, expected to be a set of transcript IDs which we should use
            for determining the reading frame around a variant. If omitted, then
            try to use all overlapping reference transcripts.

        Yields pairs of a Variant and a list of ProteinSequence objects
        """
        for (variant, overlapping_reads) in variant_and_overlapping_reads_generator:
            protein_sequences = \
                self.protein_sequences_from_variant_and_overlapping_reads(
                    variant=variant,
                    overlapping_reads=overlapping_reads,
                    transcript_id_whitelist=transcript_id_whitelist)
            yield variant, protein_sequences[:self.max_protein_sequences_per_variant]
