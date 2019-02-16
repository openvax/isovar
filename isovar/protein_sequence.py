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
ProteinSequence is a representation of a translated coding sequence,
associated with its supporting (and non-supporting but overlapping) RNA reads.
"""

from __future__ import print_function, division, absolute_import

from .translation import TranslationKey
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

    @property
    def ref_read_names(self):
        return {r.name for r in self.ref_reads}

    @property
    def alt_read_names(self):
        return {r.name for r in self.alt_reads}

    @property
    def overlapping_read_names(self):
        return {r.name for r in self.overlapping_reads}

    @property
    def ref_and_alt_read_names(self):
        return self.ref_read_names.union(self.alt_read_names)

    @property
    def other_nonref_read_names(self):
        return self.overlapping_read_names.difference(self.ref_and_alt_read_names)

    @property
    def num_ref_reads(self):
        return len(self.ref_reads)

    @property
    def num_ref_fragments(self):
        return len(self.ref_read_names)

    @property
    def num_alt_reads(self):
        return len(self.alt_reads)

    @property
    def num_alt_fragments(self):
        return len(self.alt_read_names)

    @property
    def num_overlapping_reads(self):
        return len(self.overlapping_reads)

    @property
    def num_overlapping_fragments(self):
        return len(self.overlapping_read_names)

    @property
    def num_other_nonref_reads(self):
        num_nonref = self.num_overlapping_reads - self.num_ref_reads
        return num_nonref - self.num_alt_reads

    @property
    def num_other_nonref_fragments(self):
        num_nonref = self.num_overlapping_fragments - self.num_ref_fragments
        return num_nonref - self.num_alt_fragments
