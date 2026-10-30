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

from collections.abc import Iterator

from varcode import MutationEffect, iter_candidate_effects, iter_candidate_paths
from varcode.effect_hypotheses import RealizedEffectCandidate

from .logging import get_logger

logger = get_logger(__name__)


def effect_in_coding_sequence(effect: MutationEffect | None) -> bool:
    """Keep potentially coding effects, including unresolved alternatives.

    The selected effect's own flag remains unchanged. Ordinary Varcode results
    and legacy outcome sets both retain all candidate comparisons.
    """
    if effect is None:
        return False
    return effect.modifies_coding_sequence is not False or any(
        candidate is None or candidate.modifies_coding_sequence is not False
        for candidate in iter_candidate_effects(effect))


def effect_outcomes(effect: MutationEffect | None) -> Iterator[MutationEffect | None]:
    """Yield terminal candidate effects once by identity, retaining unknowns.

    This compatibility name uses Varcode's shared traversal. For evidence or
    patient-baseline comparisons use ``iter_candidate_paths`` instead, since
    distinct provenance paths may lead to the same effect.
    """
    yield from iter_candidate_effects(effect)


def effect_protein_comparisons(
        effect: MutationEffect | None,
) -> Iterator[tuple[MutationEffect | None, str | None, str | None]]:
    """Yield each comparison's effect, baseline protein and mutant protein.

    Realized outcomes retain their own patient products, including separate
    phase/baseline hypotheses which share a terminal effect. Legacy and bare
    effects use their standard sequence accessors. Unavailable products stay
    ``None``; an empty string is a known empty protein. The original candidate
    graph, evidence and ordering are never modified.
    """
    for terminal, path in iter_candidate_paths(effect):
        candidate = path[-1] if path else None
        if terminal is not None and isinstance(candidate, RealizedEffectCandidate):
            if not candidate.outcomes:
                yield None, None, None
            for outcome in candidate.outcomes:
                yield (
                    outcome.effect,
                    outcome.baseline.to_mutant_transcript().mutant_protein_sequence,
                    outcome.mutant.to_mutant_transcript().mutant_protein_sequence)
        else:
            yield (terminal, getattr(terminal, "original_protein_sequence", None),
                   getattr(terminal, "mutant_protein_sequence", None))


def has_mutant_protein_sequence(effect: MutationEffect | None) -> bool:
    """Whether at least one outcome has a known protein (possibly empty)."""
    return any(mutant is not None for _, _, mutant in effect_protein_comparisons(effect))

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
        If True, retain effects which can change coding sequence, including
        silent changes, potential splice consequences and unresolved effects.
        Outcome sets are retained intact.

    require_mutant_protein_sequence : bool
        Drop effects without any known candidate protein sequence. Retained
        sets can still contain unresolved alternatives.

    Returns a varcode.EffectCollection object
    """
    effects = variant.effects(raise_on_error=False)
    n_total_effects = len(effects)
    logger.debug("Predicted total %d effects for variant %s", n_total_effects, variant)

    # effects filtered by allowed transcripts
    effects_filtered_by_transcript = []

    for effect in effects:
        has_transcript = hasattr(effect, 'transcript') and effect.transcript is not None
        transcript = getattr(effect, 'transcript', None)
        transcript_is_coding = (
            has_transcript and
            transcript.is_protein_coding and
            transcript.complete)
        if only_coding_transcripts and not transcript_is_coding:
            continue
        elif transcript_id_whitelist is not None and not has_transcript:
            continue
        elif transcript_id_whitelist is not None and transcript.id not in transcript_id_whitelist:
            logger.debug(
                "Skipping transcript %s for variant %s because it's not in whitelist",
                transcript.name,
                variant)
            continue
        effects_filtered_by_transcript.append(effect)

    effects = effects.clone_with_new_elements(effects_filtered_by_transcript)
    logger.debug(
        "Keeping %d/%d effects which have associated coding transcripts for %s: %s",
        len(effects),
        n_total_effects,
        variant,
        effects)

    if only_coding_effects:
        effects = effects.clone_with_new_elements([
            e for e in effects if effect_in_coding_sequence(e)
        ])
        logger.debug(
            "Keeping %d/%d effects in coding sequence for %s: %s",
            len(effects),
            n_total_effects,
            variant,
            effects)

    if require_mutant_protein_sequence:
        effects_with_mut_sequence = [
            effect
            for effect in effects
            if has_mutant_protein_sequence(effect)
        ]
        logger.debug(
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
    variant and for which it can affect the coding sequence, including
    unresolved splice outcomes and synonymous changes.

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
