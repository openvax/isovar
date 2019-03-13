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
Helper functions used to combine variant cDNA sequences collected from a
BAM file with the reading frames of annotated reference transcripts to create
candidate translations.
"""


from __future__ import print_function, division, absolute_import


from .logging import get_logger
from .reference_context import reference_contexts_for_variant
from .variant_sequence_helpers import reads_to_variant_sequences

logger = get_logger(__name__)

def translation_generator(
        variant_sequences,
        reference_contexts,
        min_transcript_prefix_length,
        max_transcript_mismatches,
        include_mismatches_after_variant,
        protein_sequence_length=None):
    """
    Given all detected VariantSequence objects for a particular variant
    and all the ReferenceContext objects for that locus, translate
    multiple protein sequences, up to the number specified by the argument
    max_protein_sequences_per_variant.

    Parameters
    ----------
    variant_sequences : list of VariantSequence objects
        Variant sequences overlapping a single original variant

    reference_contexts : list of ReferenceContext objects
        Reference sequence contexts from the same variant as the variant_sequences

    min_transcript_prefix_length : int
        Minimum number of nucleotides before the variant to test whether
        our variant sequence can use the reading frame from a reference
        transcript.

    max_transcript_mismatches : int
        Maximum number of mismatches between coding sequence before variant
        and reference transcript we're considering for determing the reading
        frame.

    include_mismatches_after_variant : bool
        If true, mismatches occurring after the variant locus will also count
        toward max_transcript_mismatches filtering.

    protein_sequence_length : int, optional
        Truncate protein to be at most this long.

    Yields a sequence of Translation objects.
    """
    for reference_context in reference_contexts:
        for variant_sequence in variant_sequences:
            translation = Translation.from_variant_sequence_and_reference_context(
                variant_sequence=variant_sequence,
                reference_context=reference_context,
                min_transcript_prefix_length=min_transcript_prefix_length,
                max_transcript_mismatches=max_transcript_mismatches,
                include_mismatches_after_variant=include_mismatches_after_variant,
                protein_sequence_length=protein_sequence_length)
            if translation is not None:
                yield translation

