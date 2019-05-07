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
        # AlleleRead objects used to construct coding sequence used
        # to create this ProteinSequence
        "reads_supporting_protein_sequence",
        # IDs of reference transcripts used to establish the reading frame for
        # this protein sequence
        "transcripts_supporting_protein_sequence",
    ]

    def __init__(
            self,
            amino_acids,
            variant_aa_interval_start,
            variant_aa_interval_end,
            ends_with_stop_codon,
            frameshift,
            translations,
            reads_supporting_protein_sequence,
            transcripts_supporting_protein_sequence):
        # fields of TranslationKey
        self.amino_acids = amino_acids
        self.variant_aa_interval_start = variant_aa_interval_start
        self.variant_aa_interval_end = variant_aa_interval_end
        self.ends_with_stop_codon = ends_with_stop_codon
        self.frameshift = frameshift

        # extra fields added by ProteinSequence
        self.translations = translations
        self.reads_supporting_protein_sequence = reads_supporting_protein_sequence
        self.transcripts_supporting_protein_sequence = transcripts_supporting_protein_sequence


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
            reads_supporting_protein_sequence,
            transcripts_supporting_protein_sequence):
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
            reads_supporting_protein_sequence=reads_supporting_protein_sequence,
            transcripts_supporting_protein_sequence=transcripts_supporting_protein_sequence)

    @property
    def num_supporting_fragments(self):
        """
        Number of unique read names used to construct the cDNA sequences from
        which this protein sequence was translated.

        Returns int
        """
        return len({r.name for r in self.reads_supporting_protein_sequence})

    @property
    def num_supporting_reads(self):
        """
        Number of reads used to construct the cDNA sequences from
        which this protein sequence was translated.

        Returns int
        """
        return len(self.reads_supporting_protein_sequence)

    @property
    def num_mismatches_before_variant(self):
        """
        Since a ProteinSequence may arise from multiple equivalent translations,
        take the minimum mismatch score from all the translations.

        Returns int
        """
        return min(t.number_mismatches_before_variant for t in self.translations)


    @property
    def num_mismatches_after_variant(self):
        """
        Since a ProteinSequence may arise from multiple equivalent translations,
        take the minimum mismatch score from all the translations.

        Returns int
        """
        return min(t.number_mismatches_after_variant for t in self.translations)

    @property
    def num_mismatches(self):
        """
        Add up the mismatches before and after the variant across all
        translations used to create this ProteinSequence.

        Returns int
        """
        return self.num_mismatches_before_variant + self.num_mismatches_after_variant

    def ascending_sort_key(self):
        """
        Sort key function used to sort protein sequences lexicographically by these criteria:
            - number of unique supporting fragments
            - number of unique supporting reads (either 1 or 2 per fragment)
            - minimum mismatch versus a supporting reference transcript before variant
            - minimum mismatch versus a supporting reference transcript after variant
            - number of supporting reference transcripts
        """
        return (
            self.num_supporting_fragments,
            self.num_supporting_reads,
            self.num_mismatches_before_variant,
            self.num_mismatches_after_variant,
            self.num_supporting_transcripts,
        )