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

from varcode import EffectCollection

from .logging import get_logger

logger = get_logger(__name__)


def predicted_effects_for_variant(
        variant,
        transcript_id_whitelist=None,
        only_coding_transcripts=False,
        drop_silent_and_noncoding=False,
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

    drop_silent_and_noncoding : bool
        If True, drop effects which aren't predicted to change the protein
        sequence.

    require_mutant_protein_sequence : bool
        Drop effects for which we can't predict what the new protein sequence
        will be.

    Returns a varcode.EffectCollection object
    """
    effects = []
    for transcript in variant.transcripts:
        if (only_coding_transcripts and not (
                transcript.complete and transcript.is_protein_coding)):
            continue
        if transcript_id_whitelist and transcript.id not in transcript_id_whitelist:
            logger.info(
                "Skipping transcript %s for variant %s because it's not in whitelist",
                transcript.name,
                variant)
            continue
        effects.append(variant.effect_on_transcript(transcript))

    effects = EffectCollection(effects)

    n_total_effects = len(effects)
    logger.info("Predicted total %d effects for variant %s" % (
        n_total_effects,
        variant))
    if drop_silent_and_noncoding:
        nonsynonymous_coding_effects = effects.drop_silent_and_noncoding()
        logger.info(
            "Keeping %d/%d effects which affect protein coding sequence for %s: %s",
            len(nonsynonymous_coding_effects),
            n_total_effects,
            variant,
            nonsynonymous_coding_effects)
        effects = nonsynonymous_coding_effects

    if require_mutant_protein_sequence:
        effects_with_mut_sequence = [
            effect
            for effect in nonsynonymous_coding_effects
            if effect.mutant_protein_sequence is not None
        ]
        logger.info(
            "Keeping %d effects with predictable AA sequences for %s: %s",
            len(effects_with_mut_sequence),
            variant,
            effects_with_mut_sequence)
        effects = effects_with_mut_sequence
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
    """
    predicted_effects = predicted_effects_for_variant(
        variant=variant,
        transcript_id_whitelist=transcript_id_whitelist,
        only_coding_transcripts=True,
        drop_silent_and_noncoding=True,
        require_mutant_protein_sequence=True)
    return [effect.transcript for effect in predicted_effects]

