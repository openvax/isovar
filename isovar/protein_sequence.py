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
        "transcript_ids_supporting_protein_sequence",
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
            transcript_ids_supporting_protein_sequence):
        # fields of TranslationKey
        self.amino_acids = amino_acids
        self.variant_aa_interval_start = variant_aa_interval_start
        self.variant_aa_interval_end = variant_aa_interval_end
        self.ends_with_stop_codon = ends_with_stop_codon
        self.frameshift = frameshift

        # extra fields added by ProteinSequence
        self.translations = translations
        self.reads_supporting_protein_sequence = reads_supporting_protein_sequence
        self.transcript_ids_supporting_protein_sequence = transcript_ids_supporting_protein_sequence

    @classmethod
    def from_translation_key(
            cls,
            translation_key,
            translations,
            reads_supporting_protein_sequence,
            transcript_ids_supporting_protein_sequence):
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
            transcript_ids_supporting_protein_sequence=transcript_ids_supporting_protein_sequence)

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
        return min(t.num_mismatches_before_variant for t in self.translations)

    @property
    def num_mismatches_after_variant(self):
        """
        Since a ProteinSequence may arise from multiple equivalent translations,
        take the minimum mismatch score from all the translations.

        Returns int
        """
        return min(t.num_mismatches_after_variant for t in self.translations)

    @property
    def num_mismatches(self):
        """
        Add up the mismatches before and after the variant across all
        translations used to create this ProteinSequence.

        Returns int
        """
        return self.num_mismatches_before_variant + self.num_mismatches_after_variant

    def transcripts(self):
        """
        Ensembl transcript IDs of all transcripts which support the reading
        frame used by Translation objects associated with this
        ProteinSequence.

        Returns list of str
        """
        transcript_set = set([])
        for t in self.translations:
            transcript_set.update(t.transcript_ids_supporting_protein_sequence)
        return sorted(transcript_set)

    def transcripts_from_protein_sequences(self, protein_sequence_limit=None):
        """
        Ensembl transcripts which support the reading frame used by protein
        sequences in this IsovarResult.

        Parameters
        ----------
        protein_sequence_limit : int or None
            If supplied then only consider the top protein sequences up to
            this number.

        Returns list of pyensembl.Transcript
        """
        genome = self.variant.genome
        transcript_ids = self.transcript_ids_used_by_protein_sequences(
            num_protein_sequences=protein_sequence_limit)
        return [
            genome.transcript_by_id(transcript_id)
            for transcript_id in transcript_ids
        ]

    def genes_from_protein_sequences(self, protein_sequence_limit=None):
        """
        Ensembl genes which support the reading frame used by protein
        sequences in this IsovarResult.

        Parameters
        ----------
        protein_sequence_limit : int or None
            If supplied then only consider the top protein sequences up to
            this number.

        Returns list of pyensembl.Gene
        """
        transcripts = self.transcripts_used_by_protein_sequences(
            protein_sequence_limit=protein_sequence_limit)
        genes = [t.gene for t in transcripts]
        return sorted(genes)

    def gene_ids_from_protein_sequences(self, protein_sequence_limit=None):
        """
        Ensembl genes IDs which support the reading frame used by protein
        sequences in this IsovarResult.

        Parameters
        ----------
        protein_sequence_limit : int or None
            If supplied then only consider the top protein sequences up to
            this number.

        Returns list of str
        """
        return [
            g.id
            for g
            in self.genes_from_protein_sequences(protein_sequence_limit=None)
        ]

    def ascending_sort_key(self):
        """
        Sort key function used to sort protein sequences lexicographically by these criteria:
            - number of unique supporting fragments
            - number of unique supporting reads (either 1 or 2 per fragment)
            - minimum mismatch versus a supporting reference transcript before variant
            - minimum mismatch versus a supporting reference transcript after variant
            - all else being equal, prefer longer sequences
        """
        return (
            self.num_supporting_fragments,
            self.num_supporting_reads,
            -self.num_mismatches_before_variant,
            -self.num_mismatches_after_variant,
            len(self.amino_acids),
        )