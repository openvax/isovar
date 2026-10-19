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

from types import SimpleNamespace

from isovar.effect_prediction import (
    predicted_effects_for_variant,
    top_varcode_effect,
    reference_coding_transcripts_for_variant
)
from varcode import Variant
from varcode.effects import Intergenic
from .common import eq_
from .genomes_for_testing import grch38
from isovar.reference_context_helpers import reference_contexts_for_variant
from isovar.isovar_result import IsovarResult
from varcode import EffectCandidate, SpliceOutcomeSet
from varcode.effects.effect_classes import (
    ExonSkipping,
    ExonicSpliceSite,
    IntronRetention,
    NormalSplicing,
    SpliceDonor,
    Substitution,
)
from varcode.effects.effect_collection import EffectCollection

# intergenic variant from error log of using Isovar 1.0.0
intergenic_variant = Variant('1', 30256419, 'G', 'T', 'GRCh37')


def test_top_effect_outside_of_gene():
    top_effect = top_varcode_effect(intergenic_variant)
    assert top_effect is not None
    assert isinstance(top_effect, Intergenic)


def test_reference_coding_transcripts_outside_of_gene():
    transcripts = \
        reference_coding_transcripts_for_variant(
            intergenic_variant)
    eq_(len(transcripts), 0)


def test_empty_transcript_whitelist_excludes_transcriptless_effects():
    assert len(predicted_effects_for_variant(intergenic_variant, transcript_id_whitelist=set())) == 0
    assert len(predicted_effects_for_variant(intergenic_variant, transcript_id_whitelist=None)) > 0


def test_empty_transcript_whitelist_does_not_broaden_reference_scope():
    variant = Variant("14", 101980529, "G", "A", ensembl=grch38)
    unrestricted = reference_contexts_for_variant(variant, context_size=30)
    assert unrestricted
    assert reference_contexts_for_variant(variant, context_size=30, transcript_id_whitelist=set()) == []
    tid = unrestricted[0].transcripts[0].id
    restricted = reference_contexts_for_variant(variant, context_size=30, transcript_id_whitelist={tid})
    assert restricted
    assert {t.id for c in restricted for t in c.transcripts} == {tid}


def _small_splice_outcome(coding=True, resolved=True):
    """Build real Varcode outcome classes over a five-residue transcript."""
    class SmallVariant:
        def effects(self, raise_on_error=False):
            return EffectCollection([outcome], sort_key=False)

    variant = SmallVariant()
    transcript = SimpleNamespace(
        id="tiny", name="tiny", gene=SimpleNamespace(name="tiny", id="tiny"),
        complete=True, is_protein_coding=True, protein_sequence="MKTAA")
    coding_effect = Substitution(variant, transcript, 2, "T", "A") if coding else None
    normal = NormalSplicing(variant, transcript, None, coding_effect)
    alternate = ExonSkipping(
        variant, transcript, None,
        affected_exon=SimpleNamespace(exon_id="tiny-exon"), in_frame=True,
        aa_alt="QQ" if resolved else None,
        aa_mutation_start_offset=2 if resolved else None,
        mutant_protein_sequence="MKQQA" if resolved else None)
    unresolved = IntronRetention(variant, transcript, None)
    outcome = SpliceOutcomeSet(
        variant, transcript,
        candidates=tuple(EffectCandidate(e) for e in (normal, alternate, unresolved)),
        disrupted_signal_class=ExonicSpliceSite if coding else SpliceDonor)
    return variant, transcript, outcome


def test_splice_outcome_stays_whole_in_coding_effect_filter():
    variant, transcript, outcome = _small_splice_outcome()
    effects = predicted_effects_for_variant(
        variant, transcript_id_whitelist={"tiny"},
        only_coding_transcripts=True, only_coding_effects=True,
        require_mutant_protein_sequence=True)
    assert list(effects) == [outcome]
    assert len(outcome.candidates) == 3
    assert reference_coding_transcripts_for_variant(variant) == [transcript]


def test_intronic_splice_outcome_is_not_in_coding_sequence():
    variant, _, outcome = _small_splice_outcome(coding=False)
    assert outcome.effect_if_splicing_unchanged is None
    assert not predicted_effects_for_variant(variant, only_coding_effects=True)


def test_result_matches_any_resolved_splice_protein():
    variant, _, outcome = _small_splice_outcome()
    protein = SimpleNamespace(
        amino_acids="MKQQA", mutation_start_idx=2, mutation_end_idx=4)
    result = IsovarResult(variant, None, outcome, [protein])
    assert outcome.mutant_protein_sequence is None
    assert result.has_mutant_protein_sequence_from_predicted_effect
    assert result.trimmed_predicted_mutant_protein_sequence is None
    assert result.num_amino_acid_mismatches_from_predicted_effect == 0
    assert result.protein_sequence_matches_predicted_mutation_effect


def test_result_does_not_claim_unresolved_splice_protein():
    variant, _, outcome = _small_splice_outcome(coding=False, resolved=False)
    protein = SimpleNamespace(
        amino_acids="MKQQA", mutation_start_idx=2, mutation_end_idx=4)
    result = IsovarResult(variant, None, outcome, [protein])
    assert not result.has_mutant_protein_sequence_from_predicted_effect
    assert result.num_amino_acid_mismatches_from_predicted_effect is None
    assert result.protein_sequence_matches_predicted_mutation_effect is None
