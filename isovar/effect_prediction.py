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

from __future__ import print_function, division, absolute_import

from varcode import EffectCollection

from .logging import get_logger

logger = get_logger(__name__)


def predicted_effects_for_variant(
        variant,
        transcript_id_whitelist=None,
        only_coding_changes=True):
    """
    For a given variant, return its set of predicted effects. Optionally
    filter to transcripts where this variant results in a non-synonymous
    change to the protein sequence.

    Parameters
    ----------
    variant : varcode.Variant

    transcript_id_whitelist : set
        Filter effect predictions to only include these transcripts

    Returns a varcode.EffectCollection object
    """

    effects = []
    for transcript in variant.transcripts:
        if only_coding_changes and not transcript.complete:
            logger.info(
                "Skipping transcript %s for variant %s because it's incomplete",
                transcript.name,
                variant)
            continue

        if transcript_id_whitelist and transcript.id not in transcript_id_whitelist:
            logger.info(
                "Skipping transcript %s for variant %s because it's not one of %d allowed",
                transcript.name,
                variant,
                len(transcript_id_whitelist))
            continue
        effects.append(variant.effect_on_transcript(transcript))

    effects = EffectCollection(effects)

    n_total_effects = len(effects)
    logger.info("Predicted total %d effects for variant %s" % (
        n_total_effects,
        variant))
    if not only_coding_changes:
        return effects
    else:
        nonsynonymous_coding_effects = effects.drop_silent_and_noncoding()
        logger.info(
            "Keeping %d/%d effects which affect protein coding sequence for %s: %s",
            len(nonsynonymous_coding_effects),
            n_total_effects,
            variant,
            nonsynonymous_coding_effects)

        usable_effects = [
            effect
            for effect in nonsynonymous_coding_effects
            if effect.mutant_protein_sequence is not None
        ]
        logger.info(
            "Keeping %d effects with predictable AA sequences for %s: %s",
            len(usable_effects),
            variant,
            usable_effects)
        return usable_effects


def reference_transcripts_for_variant(
        variant,
        transcript_id_whitelist=None,
        only_coding_changes=True):
    """
    For a given variant, find all the transcripts which overlap the
    variant and for which it has a predictable effect on the amino acid
    sequence of the protein.
    """
    predicted_effects = predicted_effects_for_variant(
        variant=variant,
        transcript_id_whitelist=transcript_id_whitelist,
        only_coding_changes=only_coding_changes)
    return [effect.transcript for effect in predicted_effects]
