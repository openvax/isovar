"""Regression coverage for #522's filtering and RNA/prediction comparisons."""
from types import SimpleNamespace

import pytest
from varcode import EffectCandidate, MultiOutcomeEffect, SpliceOutcomeSet, Variant
from varcode.effects import EffectCollection, MutationEffect

from isovar.effect_prediction import (
    effect_in_coding_sequence, effect_outcomes, has_mutant_protein_sequence,
    predicted_effects_for_variant, reference_coding_transcripts_for_variant,
    top_varcode_effect,
)
from isovar.isovar_result import IsovarResult
from isovar.read_evidence import ReadEvidence


class Outcomes(MultiOutcomeEffect):
    def __init__(self, *effects):
        super().__init__(variant=None)
        self.candidates = tuple(EffectCandidate(e, source="fixture") for e in effects)


def prediction(sequence="MACDEFG", offset=2, aa_alt="C", coding=True):
    effect = MutationEffect(variant=None)
    effect.mutant_protein_sequence = sequence
    effect.aa_mutation_start_offset = offset
    effect.aa_alt = aa_alt
    effect._original_protein_sequence = "MAYDEFG"
    effect.modifies_coding_sequence = coding
    return effect


def result(effect, protein="ACD"):
    sequence = SimpleNamespace(amino_acids=protein, mutation_start_idx=1, mutation_end_idx=2)
    return IsovarResult(variant=None, read_evidence=None, predicted_effect=effect,
                        sorted_protein_sequences=[] if protein is None else [sequence])


def test_real_scp2_splice_set_survives_filters_and_reference_selection():
    variant = Variant("1", 53513530, "A", "C", genome=75)
    scope = {"ENST00000371514"}
    effects = predicted_effects_for_variant(
        variant, transcript_id_whitelist=scope, only_coding_effects=True,
        require_mutant_protein_sequence=True)
    effect, = effects
    assert isinstance(effect, SpliceOutcomeSet)
    assert effect.effect_if_splicing_unchanged.short_description == "p.D490A"
    assert effect.mutant_protein_sequence is None
    assert len({e.mutant_protein_sequence for e in effect_outcomes(effect)
                if e.mutant_protein_sequence is not None}) >= 2
    assert len(effect.candidates) >= 3
    assert any(e.mutant_protein_sequence is None for e in effect_outcomes(effect))
    assert {t.id for t in reference_coding_transcripts_for_variant(variant, scope)} == scope
    assert isinstance(top_varcode_effect(variant, scope), SpliceOutcomeSet)
    assert reference_coding_transcripts_for_variant(variant, set()) == []
    r = IsovarResult(variant=variant, predicted_effect=effect,
                     read_evidence=ReadEvidence.from_variant_and_allele_reads(variant, []))
    record = r.to_record()
    assert record["predicted_effect_class"] == "SpliceOutcomeSet"
    assert "exon-skip" in record["predicted_effect"]
    assert record["has_mutant_protein_sequence_from_predicted_effect"] is True
    assert record["predicted_effect_mutant_protein_sequence"] is None
    assert record["protein_sequence_matches_predicted_mutation_effect"] is None


@pytest.mark.parametrize("coding, expected", [(True, True), (False, False), (None, True)])
def test_coding_filter_keeps_uncertainty(coding, expected):
    assert effect_in_coding_sequence(prediction(coding=coding)) is expected
    assert effect_in_coding_sequence(Outcomes(prediction(coding=False),
                                              prediction(coding=coding))) is expected


def test_filters_preserve_whole_nested_set_and_provenance():
    coding = prediction()
    coding.transcript = SimpleNamespace(id="coding", complete=True, is_protein_coding=True)
    effect = Outcomes(prediction(sequence=None, coding=False), Outcomes(coding))
    effect.transcript = coding.transcript
    effect.variant = Variant("1", 1, "A", "C", genome=75)
    original = EffectCollection([effect])
    variant = SimpleNamespace(effects=lambda **kwargs: original)
    filtered = predicted_effects_for_variant(
        variant, only_coding_transcripts=True, only_coding_effects=True,
        require_mutant_protein_sequence=True)
    assert list(filtered) == [effect]
    assert filtered[0] is effect
    assert filtered[0].candidates[1].source == "fixture"
    assert len(effect.candidates) == 2


def test_self_candidates_nested_sets_and_cycles_terminate():
    effect = Outcomes()
    effect.mutant_protein_sequence = "M"
    effect.candidates = (EffectCandidate(effect), EffectCandidate(Outcomes(prediction())))
    assert list(effect_outcomes(effect))[0] is effect
    assert len(list(effect_outcomes(effect))) == 2
    assert has_mutant_protein_sequence(effect)
    cycle = Outcomes(effect)
    effect.candidates = (EffectCandidate(cycle),)
    assert list(effect_outcomes(effect)) == [None]


@pytest.mark.parametrize("sequence, expected", [(None, False), ("", True), ("M", True)])
def test_known_candidate_sequence_availability(sequence, expected):
    effect = Outcomes(prediction(sequence=sequence))
    assert has_mutant_protein_sequence(effect) is expected
    assert result(effect).has_mutant_protein_sequence_from_predicted_effect is expected


def test_match_uses_second_candidate_and_its_own_offset():
    effect = Outcomes(prediction(sequence="MAYDEFG"),
                      prediction(sequence="MMMACDEFG", offset=4),
                      prediction(sequence=None))
    r = result(effect)
    assert r.predicted_effect is effect
    assert r.trimmed_predicted_mutant_protein_sequence is None
    assert r.num_amino_acid_mismatches_from_predicted_effect == 0
    assert r.protein_sequence_matches_predicted_mutation_effect is True


@pytest.mark.parametrize("include_unknown, expected", [(False, False), (True, None)])
def test_nonmatch_is_not_negative_when_an_outcome_is_unresolved(include_unknown, expected):
    effects = [prediction(sequence="MAYDEFG"), prediction(sequence="MAWDEFG")]
    if include_unknown:
        effects.append(prediction(sequence=None))
    r = result(Outcomes(*effects))
    assert r.num_amino_acid_mismatches_from_predicted_effect == 1
    assert r.protein_sequence_matches_predicted_mutation_effect is expected


@pytest.mark.parametrize("effect", [None, Outcomes(), Outcomes(prediction(sequence=None)),
                                  Outcomes(prediction(offset=None))])
def test_no_comparable_candidate_is_unknown(effect):
    r = result(effect)
    assert r.trimmed_predicted_mutant_protein_sequence is None
    assert r.num_amino_acid_mismatches_from_predicted_effect is None
    assert r.protein_sequence_matches_predicted_mutation_effect is None


def test_scalar_windows_require_agreement_and_keep_single_effect_behavior():
    for effect in [prediction(), Outcomes(prediction(), prediction())]:
        r = result(effect)
        assert r.trimmed_predicted_mutant_protein_sequence == "ACD"
        assert r.trimmed_reference_protein_sequence == "AYD"
        assert r.protein_sequence_matches_predicted_mutation_effect is True
    assert result(Outcomes(prediction(), prediction(offset=3))).trimmed_reference_protein_sequence is None
    r = result(prediction(), protein=None)
    assert r.trimmed_predicted_mutant_protein_sequence is None
    assert r.protein_sequence_matches_predicted_mutation_effect is None
