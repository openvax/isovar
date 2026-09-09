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
