# Copyright (c) 2016. Mount Sinai School of Medicine
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
import logging

from varcode import EffectCollection


logger = logging.getLogger(__name__)


def predicted_coding_effects_with_mutant_sequence(
        variant,
        transcript_id_whitelist=None):
    """
    For a given variant, return the set of predicted mutation effects
    on transcripts where this variant results in a predictable non-synonymous
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
        if not transcript.complete:
            logger.info(
                "Skipping transcript %s for variant %s because it's incomplete",
                    transcript,
                    variant)
            continue

        if transcript_id_whitelist and transcript.id not in transcript_id_whitelist:
            logger.info(
                "Skipping transcript %s for variant %s because it's not one of %d allowed",
                    transcript,
                    variant,
                    len(transcript_id_whitelist))
            continue
        effects.append(variant.effect_on_transcript(transcript))

    effects = EffectCollection(effects)

    n_total_effects = len(effects)
    logger.info("Predicted %d effects for variant %s" % (
        n_total_effects,
        variant))

    nonsynonymous_coding_effects = effects.drop_silent_and_noncoding()
    logger.info(
        "Keeping %d/%d non-synonymous coding effects for %s" % (
            len(nonsynonymous_coding_effects),
            n_total_effects,
            variant))

    usable_effects = [
        effect
        for effect in nonsynonymous_coding_effects
        if effect.mutant_protein_sequence is not None
    ]
    logger.info(
        "Keeping %d/%d effects with predictable AA sequences for %s",
            len(usable_effects),
            len(nonsynonymous_coding_effects),
            variant)
    return usable_effects

def reference_transcripts_for_variant(variant, transcript_id_whitelist=None):
    """
    For a given variant, find all the transcripts which overlap the
    variant and for which it has a predictable effect on the amino acid
    sequence of the protein.
    """
    predicted_effects = predicted_coding_effects_with_mutant_sequence(
        variant=variant,
        transcript_id_whitelist=transcript_id_whitelist)
    return [effect.transcript for effect in predicted_effects]
