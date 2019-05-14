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

from __future__ import print_function, division, absolute_import
from collections import OrderedDict, defaultdict

from .effect_prediction import reference_coding_transcripts_for_variant
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


def reference_contexts_for_variant(
        variant,
        context_size,
        transcript_id_whitelist=None):
    """
    variant : varcode.Variant

    context_size : int
        Max of nucleotides to include to the left and right of the variant
        in the context sequence.

    transcript_id_whitelist : set, optional
        If given, then only consider transcripts whose IDs are in this set.

    Returns list of ReferenceContext objects, sorted by maximum length of
    coding sequence of any supporting transcripts.
    """
    overlapping_transcripts = reference_coding_transcripts_for_variant(
        variant=variant,
        transcript_id_whitelist=transcript_id_whitelist)

    # dictionary mapping SequenceKeyWithReadingFrame keys to list of
    # transcript objects
    sequence_groups = defaultdict(list)

    for transcript in overlapping_transcripts:
        reference_coding_sequence_key = \
            ReferenceCodingSequenceKey.from_variant_and_transcript(
                variant=variant,
                transcript=transcript,
                context_size=context_size)
        if reference_coding_sequence_key is not None:
            sequence_groups[reference_coding_sequence_key].append(transcript)

    reference_contexts = [
        ReferenceContext.from_reference_coding_sequence_key(
            key, variant, matching_transcripts)
        for (key, matching_transcripts) in sequence_groups.items()
    ]
    reference_contexts.sort(
        key=ReferenceContext.sort_key_decreasing_max_length_transcript_cds)
    return reference_contexts


def reference_contexts_for_variants(
        variants,
        context_size,
        transcript_id_whitelist=None):
    """
    Extract a set of reference contexts for each variant in the collection.

    Parameters
    ----------
    variants : varcode.VariantCollection

    context_size : int
        Max of nucleotides to include to the left and right of the variant
        in the context sequence.

    transcript_id_whitelist : set, optional
        If given, then only consider transcripts whose IDs are in this set.

    Returns a dictionary from variants to lists of ReferenceContext objects,
    sorted by max coding sequence length of any transcript.
    """
    result = OrderedDict()
    for variant in variants:
        result[variant] = reference_contexts_for_variant(
            variant=variant,
            context_size=context_size,
            transcript_id_whitelist=transcript_id_whitelist)
    return result
