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

from dataclasses import dataclass

import pytest

from varcode import (
    MutantTranscriptSource,
    ReadPhaseResolver,
    ReadPhasingSource,
    apply_phase_resolver_to_effects,
    load_vcf,
)

from isovar import IsovarMutantTranscript, IsovarReadPhasing, run_isovar

from .common import eq_
from .testing_helpers import data_path


class IsovarReadEvidence(IsovarReadPhasing, IsovarMutantTranscript):
    pass


@pytest.fixture(scope="module")
def b16_results():
    return run_isovar(
        variants=data_path("data/b16.f10/b16.vcf"),
        alignment_file=data_path("data/b16.f10/b16.combined.sorted.bam"))


@pytest.fixture(scope="module")
def result_and_transcript(b16_results):
    result = next(
        isovar_result
        for isovar_result in b16_results
        if isovar_result.has_mutant_protein_sequence_from_rna
    )
    transcript = result.top_protein_sequence.transcripts[0]
    return b16_results, result, transcript


def test_isovar_mutant_transcript_satisfies_runtime_protocol(result_and_transcript):
    results, _, _ = result_and_transcript
    source = IsovarMutantTranscript(results)
    assert isinstance(source, MutantTranscriptSource)


def test_composed_isovar_source_satisfies_read_phase_resolver_protocols(
        result_and_transcript):
    results, _, _ = result_and_transcript
    source = IsovarReadEvidence(results)
    assert isinstance(source, ReadPhasingSource)
    assert isinstance(source, MutantTranscriptSource)


def test_composed_isovar_source_initializes_both_provider_indexes(
        result_and_transcript):
    results, result, transcript = result_and_transcript
    source = IsovarReadEvidence(results)

    assert source.has_evidence(result.variant)
    assert source.mutant_transcript(result.variant, transcript) is not None


def test_isovar_mutant_transcript_exposes_observed_transcript(result_and_transcript):
    results, result, transcript = result_and_transcript
    source = IsovarMutantTranscript(results)

    mutant_transcript = source.mutant_transcript(result.variant, transcript)
    assert mutant_transcript is not None
    eq_(mutant_transcript.reference_transcript.id, transcript.id)
    assert mutant_transcript.cdna_sequence in result.top_protein_sequence.cdna_sequences
    eq_(mutant_transcript.mutant_protein_sequence, result.top_protein_sequence.amino_acids)
    eq_(mutant_transcript.annotator_name, "isovar")
    assert len(mutant_transcript.edits) >= 1
    assert any(
        edit.source_variant is not None and
        edit.source_variant == result.variant
        for edit in mutant_transcript.edits
    )


def test_read_phase_resolver_applies_isovar_mutant_transcript(result_and_transcript):
    results, result, transcript = result_and_transcript
    source = IsovarReadEvidence(results)
    resolver = ReadPhaseResolver(source)

    @dataclass
    class DummyEffect(object):
        variant: object
        transcript: object
        mutant_transcript: object = None

    effect = DummyEffect(variant=result.variant, transcript=transcript)
    effects = [effect]

    resolved_effects = apply_phase_resolver_to_effects(effects, resolver)
    eq_(resolved_effects, effects)
    assert effect.mutant_transcript is not None
    eq_(effect.mutant_transcript.reference_transcript.id, transcript.id)
    eq_(
        resolver.phased_partners(result.variant, transcript),
        source.partners_in_cis(result.variant))


def test_b16_effects_accept_composed_isovar_source(result_and_transcript):
    results, result, transcript = result_and_transcript
    source = IsovarReadEvidence(results)
    resolver = ReadPhaseResolver(source)
    tumor = load_vcf(data_path("data/b16.f10/b16.vcf"))

    effects = tumor.effects(phase_resolver=resolver)
    matching_effects = [
        effect
        for effect in effects
        if (
            getattr(effect, "variant", None) == result.variant and
            getattr(getattr(effect, "transcript", None), "id", None) == transcript.id
        )
    ]
    assert matching_effects
    mutant_transcripts = [
        effect.mutant_transcript
        for effect in matching_effects
        if getattr(effect, "mutant_transcript", None) is not None
    ]
    assert mutant_transcripts
    eq_(mutant_transcripts[0].reference_transcript.id, transcript.id)
    eq_(
        mutant_transcripts[0].mutant_protein_sequence,
        result.top_protein_sequence.amino_acids)
