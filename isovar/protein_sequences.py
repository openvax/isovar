# Copyright (c) 2016-2018. Mount Sinai School of Medicine
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

from .common import groupby
from .default_parameters import (
    MIN_TRANSCRIPT_PREFIX_LENGTH,
    MAX_REFERENCE_TRANSCRIPT_MISMATCHES,
    INCLUDE_MISMATCHES_AFTER_VARIANT,
    PROTEIN_SEQUENCE_LENGTH,
    MAX_PROTEIN_SEQUENCES_PER_VARIANT,
    MIN_ALT_RNA_READS,
    MIN_VARIANT_SEQUENCE_COVERAGE,
    VARIANT_SEQUENCE_ASSEMBLY
)
from .dataframe_builder import dataframe_from_generator
from .translation import translate_variant_reads, Translation, TranslationKey
from .read_helpers import group_reads_by_allele
from .variant_helpers import trim_variant
from .logging import get_logger

logger = get_logger(__name__)


class ProteinSequence(TranslationKey):
    """
    Translated amino acid sequence aggregated across possibly multiple
    VariantSequence and ReferenceContext objects (e.g. imagine two distinct
    sequences which contain synonymous codons).

    This is the final result of the isovar variant->expressed peptide pipeline.
    """
    __slots__ = [
        # list of all the Translation objects which support this distinct
        # amino acid sequence
        "translations",
        # number of reads overlapping the variant locus supporting any allele,
        # including the reference, alt, or anything else
        "overlapping_reads",
        # number of reads overlapping this locus which support the reference
        # allele
        "ref_reads",
        # total number of reads at the locus which contained the variant
        # nucleotides, even if they supported other phased sequences
        "alt_reads",
        # number of unique read names from all the VariantSequence objects
        # from each translation
        "alt_reads_supporting_protein_sequence",
        # IDs of transcripts overlapping the variant locus
        "transcripts_overlapping_variant",
        # IDs of reference transcripts used to establish the reading frame for
        # this protein sequence
        "transcripts_supporting_protein_sequence",
        # name of gene of the reference transcripts used in Translation
        # objects
        "gene",
    ]

    def __init__(
            self,
            amino_acids,
            variant_aa_interval_start,
            variant_aa_interval_end,
            ends_with_stop_codon,
            frameshift,
            translations,
            overlapping_reads,
            ref_reads,
            alt_reads,
            alt_reads_supporting_protein_sequence,
            transcripts_overlapping_variant,
            transcripts_supporting_protein_sequence,
            gene):
        self.amino_acids = amino_acids
        self.variant_aa_interval_start = variant_aa_interval_start
        self.variant_aa_interval_end = variant_aa_interval_end
        self.ends_with_stop_codon = ends_with_stop_codon
        self.frameshift = frameshift
        self.translations = translations
        self.overlapping_reads = overlapping_reads
        self.ref_reads = ref_reads
        self.alt_reads = alt_reads
        self.alt_reads_supporting_protein_sequence = (
            alt_reads_supporting_protein_sequence)
        self.transcripts_overlapping_variant = transcripts_overlapping_variant
        self.transcripts_supporting_protein_sequence = (
            transcripts_supporting_protein_sequence)
        self.gene = gene

    @classmethod
    def _summarize_translations(cls, translations):
        """
        Summarize a collection of Translation objects into three values:
            1) List of unique reads supporting underlying variant sequences
            2) Set of unique transcript names for establishing reading frames of the
               translations.
            3) Set of unique gene names for all transcripts used by translations.
        """
        read_name_to_reads = {}
        gene_names = set([])
        transcript_ids = set([])
        for translation in translations:
            for read in translation.reads:
                read_name_to_reads[read.name] = read
            for transcript in translation.reference_context.transcripts:
                transcript_ids.add(transcript.id)
                gene_names.add(transcript.gene.name)
        unique_reads = list(read_name_to_reads.values())
        return unique_reads, transcript_ids, gene_names

    @classmethod
    def from_translation_key(
            cls,
            translation_key,
            translations,
            overlapping_reads,
            ref_reads,
            alt_reads,
            alt_reads_supporting_protein_sequence,
            transcripts_overlapping_variant,
            transcripts_supporting_protein_sequence,
            gene):
        """
        Create a ProteinSequence object from a TranslationKey, along with
        all the extra fields a ProteinSequence requires.
        """
        return cls(
            amino_acids=translation_key.amino_acids,
            variant_aa_interval_start=translation_key.variant_aa_interval_start,
            variant_aa_interval_end=translation_key.variant_aa_interval_end,
            ends_with_stop_codon=translation_key.ends_with_stop_codon,
            frameshift=translation_key.frameshift,
            translations=translations,
            overlapping_reads=overlapping_reads,
            ref_reads=ref_reads,
            alt_reads=alt_reads,
            alt_reads_supporting_protein_sequence=(
                alt_reads_supporting_protein_sequence),
            transcripts_overlapping_variant=transcripts_overlapping_variant,
            transcripts_supporting_protein_sequence=(
                transcripts_supporting_protein_sequence),
            gene=gene)

    def ascending_sort_key(self):
        """
        Sort protein sequences lexicographically by three criteria:
            - number of unique supporting reads
            - minimum mismatch versus a supporting reference transcript before variant
            - minimum mismatch versus a supporting reference transcript after variant
            - number of supporting reference transcripts

        TODO: Add sort criterion:
            - min number of reads covering each nucleotide of
              the protein sequence >= 2
        """
        return (
            len(self.alt_reads_supporting_protein_sequence),
            min(t.number_mismatches_before_variant for t in self.translations),
            min(t.number_mismatches_after_variant for t in self.translations),
            len(self.transcripts_supporting_protein_sequence)
        )

def sort_protein_sequences(protein_sequences):
    """
    Sort protein sequences in decreasing order of priority
    """
    return list(
        sorted(
            protein_sequences,
            key=ProteinSequence.ascending_sort_key,
            reverse=True))

def reads_generator_to_protein_sequences_generator(
        variant_and_overlapping_reads_generator,
        transcript_id_whitelist=None,
        protein_sequence_length=PROTEIN_SEQUENCE_LENGTH,
        min_alt_rna_reads=MIN_ALT_RNA_READS,
        min_variant_sequence_coverage=MIN_VARIANT_SEQUENCE_COVERAGE,
        min_transcript_prefix_length=MIN_TRANSCRIPT_PREFIX_LENGTH,
        max_transcript_mismatches=MAX_REFERENCE_TRANSCRIPT_MISMATCHES,
        include_mismatches_after_variant=INCLUDE_MISMATCHES_AFTER_VARIANT,
        max_protein_sequences_per_variant=MAX_PROTEIN_SEQUENCES_PER_VARIANT,
        variant_sequence_assembly=VARIANT_SEQUENCE_ASSEMBLY):
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

    protein_sequence_length : int
        Try to translate protein sequences of this length, though sometimes
        we'll have to return something shorter (depending on the RNAseq data,
        and presence of stop codons).

    min_alt_rna_reads : int
        Drop variant sequences at loci with fewer than this number of reads
        supporting the alt allele.

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

    variant_cdna_sequence_assembly : bool
        If True, then assemble variant cDNA sequences based on overlap of
        RNA reads. If False, then variant cDNA sequences must be fully spanned
        and contained within RNA reads.

    Yields pairs of a Variant and a list of ProteinSequence objects
    """

    for (variant, overlapping_reads) in variant_and_overlapping_reads_generator:
        overlapping_transcript_ids = [
            t.id
            for t in variant.transcripts
            if t.is_protein_coding
        ]
        _, ref, alt = trim_variant(variant)
        overlapping_reads = list(overlapping_reads)
        reads_grouped_by_allele = group_reads_by_allele(overlapping_reads)

        ref_reads = reads_grouped_by_allele.get(ref, [])
        alt_reads = reads_grouped_by_allele.get(alt, [])

        translations = translate_variant_reads(
            variant=variant,
            variant_reads=alt_reads,
            transcript_id_whitelist=transcript_id_whitelist,
            protein_sequence_length=protein_sequence_length,
            min_alt_rna_reads=min_alt_rna_reads,
            min_variant_sequence_coverage=min_variant_sequence_coverage,
            min_transcript_prefix_length=min_transcript_prefix_length,
            max_transcript_mismatches=max_transcript_mismatches,
            include_mismatches_after_variant=include_mismatches_after_variant,
            variant_sequence_assembly=variant_sequence_assembly)

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
                alt_reads=alt_reads,
                ref_reads=ref_reads,
                alt_reads_supporting_protein_sequence=alt_reads_supporting_protein_sequence,
                transcripts_supporting_protein_sequence=group_transcript_ids,
                transcripts_overlapping_variant=overlapping_transcript_ids,
                gene=list(group_gene_names))
            logger.info("%s: protein sequence = %s" % (key, protein_sequence.amino_acids))
            protein_sequences.append(protein_sequence)

        # sort protein sequences before returning the top results
        protein_sequences = sort_protein_sequences(protein_sequences)

        yield variant, protein_sequences[:max_protein_sequences_per_variant]


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
