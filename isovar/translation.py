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
This module combines variant cDNA sequences collected from a BAM file with
the reading frames of annotated reference transcripts to create candidate
translations.
"""


from __future__ import print_function, division, absolute_import

from .logging import get_logger
from .translation_key import TranslationKey


logger = get_logger(__name__)


class Translation(TranslationKey):
    """
    Translated amino acid sequence of a VariantORF for a particular
    ReferenceContext and VariantSequence.
    """
    __slots__ = [
        "untrimmed_variant_sequence",
        "reference_context",
        "variant_orf"
    ]

    def __init__(
            self,
            amino_acids,
            contains_mutation,
            mutation_start_idx,
            mutation_end_idx,
            ends_with_stop_codon,
            frameshift,
            untrimmed_variant_sequence,
            reference_context,
            variant_orf):
        """
        Parameters
        ----------
        amino_acids : str
            Amino acid sequence

        contains_mutation : bool
            Does the amino acid sequence contain a mutation?

        mutation_start_idx : int
            Start of half-open interval for variant amino acids
            in the translated sequence

        mutation_end_idx : int
            End of half-open interval for variant amino acids
            in the translated sequence

        ends_with_stop_codon : bool
            Did the amino acid sequence end due to a stop codon or did we
            just run out of sequence context around the variant?

        frameshift : bool
            Was the variant a frameshift relative to the reference sequence?

        untrimmed_variant_sequence : VariantSequence

        reference_context : ReferenceContext

        variant_orf : VariantORF
        """

        # TODO:
        #  get rid of untrimmed_variant_sequence by making
        #  VariantORF keep track of its inputs

        self.amino_acids = amino_acids
        self.contains_mutation = contains_mutation
        self.mutation_start_idx = mutation_start_idx
        self.mutation_end_idx = mutation_end_idx
        self.ends_with_stop_codon = ends_with_stop_codon
        self.frameshift = frameshift
        # this variant sequence might differ from the one
        # in variant_orf due to trimming required to match the reference
        self.untrimmed_variant_sequence = untrimmed_variant_sequence
        self.reference_context = reference_context
        self.variant_orf = variant_orf

    @property
    def reads(self):
        """
        RNA reads which were used to construct the coding sequence
        from which we translated these amino acids.
        """
        return self.untrimmed_variant_sequence.reads

    @property
    def reference_cdna_sequence_before_variant(self):
        """

        Returns str
        """
        return self.variant_orf.reference_cdna_sequence_before_variant


    @property
    def num_mismatches_before_variant(self):
        """
        Number of nucleotides in the variant cDNA sequence which
        don't match the ReferenceContext transcript sequence at
        positions before the variant locus.

        Returns int
        """
        return self.variant_orf.num_mismatches_before_variant

    @property
    def num_mismatches_after_variant(self):
        """
        Number of nucleotides in the variant cDNA sequence which
        don't match the ReferenceContext transcript sequence at
        positions after the variant locus.

        Returns int
        """
        return self.variant_orf.num_mismatches_after_variant

    @property
    def cdna_sequence(self):
        """
        cDNA sequence assembled from variant supporting reads

        Returns str
        """
        return self.variant_orf.cdna_sequence

    @property
    def offset_to_first_complete_codon(self):
        """
        Offset to first complete codon in the cDNA sequence

        Returns int in {0, 1, 2}
        """
        return self.variant_orf.offset_to_first_complete_codon

    @property
    def variant_cdna_interval_start(self):
        """
        Interbase start coordinate of variant interval in the cDNA sequence
        """
        return self.variant_orf.variant_cdna_interval_start

    @property
    def variant_cdna_interval_end(self):
        """
        Interbase end coordinate of variant interval in the cDNA sequence
        """
        return self.variant_orf.variant_cdna_interval_end

    def as_translation_key(self):
        """
        Project Translation object or any other derived class into just a
        TranslationKey, which has fewer fields and can be used as a
        dictionary key.
        """
        return TranslationKey(**{
            name: getattr(self, name)
            for name in TranslationKey._fields})


