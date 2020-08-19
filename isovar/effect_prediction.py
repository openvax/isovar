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

from varcode.effects.effect_classes import ExonicSpliceSite, CodingMutation

from .logging import get_logger

logger = get_logger(__name__)


def effect_in_coding_sequence(effect):
    if type(effect) is ExonicSpliceSite:
        return effect_in_coding_sequence(effect.alternate_effect)
    else:
        return isinstance(effect, CodingMutation)

def predicted_effects_for_variant(
        variant,
        transcript_id_whitelist=None,
        only_coding_transcripts=False,
        only_coding_effects=False,
        require_mutant_protein_sequence=False):
    """
    For a given variant, return its set of predicted effects. Optionally
    filter to transcripts where this variant results in a non-synonymous
    change to the protein sequence.

    Parameters
    ----------
    variant : varcode.Variant

    transcript_id_whitelist : set
        Filter effect predictions to only include these transcripts

    only_coding_transcripts : bool
        If True, then only return effects on protein coding transcripts

    only_coding_effects : bool
        If True, drop effects which aren't in the coding sequence, even
        if they are silent.

    require_mutant_protein_sequence : bool
        Drop effects for which we can't predict what the new protein sequence
        will be.

    Returns a varcode.EffectCollection object
    """
    effects = variant.effects(raise_on_error=False)
    n_total_effects = len(effects)
    logger.info("Predicted total %d effects for variant %s" % (
        n_total_effects,
        variant))

    # effects filtered by allowed transcripts
    effects_filtered_by_transcript = []

    for effect in effects:
        has_transcript = hasattr(effect, 'transcript') and effect.transcript is not None
        transcript = getattr(effect, 'transcript', None)
        transcript_is_coding = (
            has_transcript and
            transcript.complete and
            transcript.is_protein_coding)
        if only_coding_transcripts and not transcript_is_coding:
            continue
        elif transcript_id_whitelist and not has_transcript:
            continue
        elif transcript_id_whitelist and transcript.id not in transcript_id_whitelist:
            logger.info(
                "Skipping transcript %s for variant %s because it's not in whitelist",
                transcript.name,
                variant)
            continue
        effects_filtered_by_transcript.append(effect)

    effects = effects.clone_with_new_elements(effects_filtered_by_transcript)
    logger.info(
        "Keeping %d/%d effects which have associated coding transcripts for %s: %s",
        len(effects),
        n_total_effects,
        variant,
        effects)

    if only_coding_effects:
        effects = effects.clone_with_new_elements([
            e for e in effects if effect_in_coding_sequence(e)
        ])
        logger.info(
            "Keeping %d/%d effects in coding sequence for %s: %s",
            len(effects),
            n_total_effects,
            variant,
            effects)

    if require_mutant_protein_sequence:
        effects_with_mut_sequence = [
            effect
            for effect in effects
            if effect.mutant_protein_sequence is not None
        ]
        logger.info(
            "Keeping %d effects with predictable AA sequences for %s: %s",
            len(effects_with_mut_sequence),
            variant,
            effects_with_mut_sequence)
        effects = effects.clone_with_new_elements(effects_with_mut_sequence)
    return effects


def top_varcode_effect(variant, transcript_id_whitelist=None):
    """
    Find the best predicted effect for the given variant. If we have a
    transcript whitelist (based on filtering bulk expression) then use
    it to eliminate some of the effect predictions.
    Returns subclass of varcode.MutationEffect
    """
    effects = predicted_effects_for_variant(
        variant,
        transcript_id_whitelist=transcript_id_whitelist)
    if len(effects) == 0 and transcript_id_whitelist is not None:
        # if everything got filtered due to the transcript whitelist,
        # we still need to return some kind of "top" effect so look
        # at those which got filtered out by expression
        effects = predicted_effects_for_variant(
            variant,
            transcript_id_whitelist=None)
    if len(effects) == 0:
        raise ValueError(
            "Could not determine top effect prediction for %s" % variant)
    return effects.top_priority_effect()


def reference_coding_transcripts_for_variant(
        variant,
        transcript_id_whitelist=None):
    """
    For a given variant, find all the transcripts which overlap the
    variant and for which it has a predictable effect on the amino acid
    sequence of the protein.

    Parameters
    ----------
    variant : varcode.Variant

    transcript_id_whitelist : set or None
    """
    predicted_effects = predicted_effects_for_variant(
        variant=variant,
        transcript_id_whitelist=transcript_id_whitelist,
        only_coding_transcripts=True,
        only_coding_effects=True,
        require_mutant_protein_sequence=False)
    return [effect.transcript for effect in predicted_effects]

