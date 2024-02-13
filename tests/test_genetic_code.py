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

from isovar.genetic_code import translate_cdna
from .common import eq_
from .genomes_for_testing import grch38


def test_translate_cdna_no_stop_codon():
    eq_(translate_cdna("ATGATG", first_codon_is_start=False), ("MM", False))


def test_translate_cdna_stop_codon():
    eq_(translate_cdna("ATGATGTAG", first_codon_is_start=False), ("MM", True))


def test_translate_cdna_alternate_CTG_start():
    eq_(translate_cdna("CTGCTG", first_codon_is_start=True), ("ML", False))


def test_translate_cdna_CTG_after_start():
    eq_(translate_cdna("CTGCTG", first_codon_is_start=False), ("LL", False))


def test_TP53_translation_from_cdna():
    tp53_001 = grch38.transcripts_by_name("TP53-001")[0]
    cdna = tp53_001.coding_sequence
    amino_acids, ends_with_stop_codon = translate_cdna(cdna, first_codon_is_start=True)
    assert ends_with_stop_codon
    eq_(amino_acids, tp53_001.protein_sequence)


def test_mitochondrial_MTND5_translation_from_cdna():
    mtnd5_001 = grch38.transcripts_by_name("MT-ND5-201")[0]
    cdna = mtnd5_001.coding_sequence
    amino_acids, ends_with_stop_codon = translate_cdna(
        cdna,
        first_codon_is_start=True,
        mitochondrial=True)
    assert ends_with_stop_codon
    eq_(amino_acids, mtnd5_001.protein_sequence)
