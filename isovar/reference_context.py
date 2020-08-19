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

from .reference_coding_sequence_key import ReferenceCodingSequenceKey
from .logging import get_logger

logger = get_logger(__name__)


##########################
#
# ReferenceContext
# ----------------
#
# Includes all the fields of ReferenceCodingSequenceKey in addition to which
# variant we're examining and all transcripts overlapping that variant
# which produced this particular sequence context and reading frame.
#
##########################

class ReferenceContext(ReferenceCodingSequenceKey):
    """
    Representation of the sequence preceding a variant and its associated
    reading frame.
    """

    # additional fields on top of slots for ReferenceCodingSequenceKey
    __slots__ = ["variant", "transcripts"]

    def __init__(
            self,
            strand,
            sequence_before_variant_locus,
            sequence_at_variant_locus,
            sequence_after_variant_locus,
            offset_to_first_complete_codon,
            contains_start_codon,
            overlaps_start_codon,
            contains_five_prime_utr,
            amino_acids_before_variant,
            variant,
            transcripts):
        ReferenceCodingSequenceKey.__init__(
            self,
            strand=strand,
            sequence_before_variant_locus=sequence_before_variant_locus,
            sequence_at_variant_locus=sequence_at_variant_locus,
            sequence_after_variant_locus=sequence_after_variant_locus,
            offset_to_first_complete_codon=offset_to_first_complete_codon,
            contains_start_codon=contains_start_codon,
            overlaps_start_codon=overlaps_start_codon,
            contains_five_prime_utr=contains_five_prime_utr,
            amino_acids_before_variant=amino_acids_before_variant)
        self.variant = variant
        self.transcripts = tuple(transcripts)

    @classmethod
    def from_reference_coding_sequence_key(cls, key, variant, transcripts):
        """
        Construct a ReferenceContext object from a ReferenceSequenceKey, variant,
        and a set of transcript.

        Parameters
        ----------
        key : ReferenceSequenceKey
        variant : varcode.Variant
        transcripts : list of pyensembl.Transcript

        Returns ReferenceContext
        """
        return ReferenceContext(
            strand=key.strand,
            sequence_before_variant_locus=key.sequence_before_variant_locus,
            sequence_at_variant_locus=key.sequence_at_variant_locus,
            sequence_after_variant_locus=key.sequence_after_variant_locus,
            offset_to_first_complete_codon=key.offset_to_first_complete_codon,
            contains_start_codon=key.contains_start_codon,
            overlaps_start_codon=key.overlaps_start_codon,
            contains_five_prime_utr=key.contains_five_prime_utr,
            amino_acids_before_variant=key.amino_acids_before_variant,
            variant=variant,
            transcripts=transcripts)

    def sort_key_decreasing_max_length_transcript_cds(self):
        """
        Used to sort a sequence of ReferenceContext objects by the longest CDS
        in each context's list of transcripts.
        """
        return -max(len(t.coding_sequence) for t in self.transcripts)

    @property
    def mitochondrial(self):
        """
        Is this a reference context for a variant in the mitochondrial
        genome?
        """
        return self.variant.contig.lower() in {"chrm", "m", "chrmt", "mt"}

