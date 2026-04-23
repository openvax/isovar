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

from varcode.phasing import (
    IsovarAssemblyProvider,
    IsovarPhaseResolver,
    apply_phase_resolver_to_effects,
)

from isovar import VarcodeAdapter, run_isovar

from .common import eq_
from .testing_helpers import data_path


@pytest.fixture(scope="module")
def result_and_transcript():
    results = run_isovar(
        variants=data_path("data/b16.f10/b16.vcf"),
        alignment_file=data_path("data/b16.f10/b16.combined.sorted.bam"))
    result = next(
        isovar_result
        for isovar_result in results
        if isovar_result.has_mutant_protein_sequence_from_rna
    )
    transcript = result.top_protein_sequence.transcripts[0]
    return results, result, transcript


def test_varcode_adapter_satisfies_runtime_protocol(result_and_transcript):
    results, _, _ = result_and_transcript
    provider = VarcodeAdapter(results)
    assert isinstance(provider, IsovarAssemblyProvider)


def test_varcode_adapter_exposes_contig_and_mutant_transcript(result_and_transcript):
    results, result, transcript = result_and_transcript
    provider = VarcodeAdapter(results)

    assert provider.has_contig(result.variant, transcript)

    mutant_transcript = provider.mutant_transcript(result.variant, transcript)
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


def test_isovar_phase_resolver_applies_contig_mutant_transcript(result_and_transcript):
    results, result, transcript = result_and_transcript
    provider = VarcodeAdapter(results)
    resolver = IsovarPhaseResolver(provider)

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
    eq_(resolver.phased_partners(result.variant, transcript), provider.variants_in_contig(result.variant, transcript))
