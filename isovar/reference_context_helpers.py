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

from collections import defaultdict


from .effect_prediction import reference_coding_transcripts_for_variant
from .reference_context import ReferenceContext
from .reference_coding_sequence_key import ReferenceCodingSequenceKey


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


def reference_contexts_generator(
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

    Generate a series of (Variant, [ReferenceContext]) pairs, where the
    to list of ReferenceContext objects for each variant is sorted by
    max coding sequence length of any transcript.
    """
    for variant in variants:
        reference_contexts = reference_contexts_for_variant(
            variant=variant,
            context_size=context_size,
            transcript_id_whitelist=transcript_id_whitelist)
        yield variant, reference_contexts