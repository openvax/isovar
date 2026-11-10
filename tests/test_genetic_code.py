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

from itertools import product

import pytest

from isovar.genetic_code import standard_genetic_code, translate_cdna
from .common import eq_
from .genomes_for_testing import grch38


def test_translate_cdna_no_stop_codon():
    eq_(translate_cdna("ATGATG", first_codon_is_start=False), ("MM", False))


def test_translate_cdna_stop_codon():
    eq_(translate_cdna("ATGATGTAG", first_codon_is_start=False, translation_table=1), ("MM", True))


def test_translate_cdna_alternate_CTG_start():
    eq_(translate_cdna("CTGCTG", first_codon_is_start=True, translation_table=1), ("ML", False))


def test_translate_cdna_CTG_after_start():
    eq_(translate_cdna("CTGCTG", first_codon_is_start=False, translation_table=1), ("LL", False))


def test_TP53_translation_from_cdna():
    tp53_001 = grch38.transcripts_by_name("TP53-001")[0]
    cdna = tp53_001.coding_sequence
    amino_acids, ends_with_stop_codon = translate_cdna(cdna, first_codon_is_start=True, translation_table=1)
    assert ends_with_stop_codon
    eq_(amino_acids, tp53_001.protein_sequence)


def test_mitochondrial_MTND5_translation_from_cdna():
    mtnd5_001 = grch38.transcripts_by_name("MT-ND5-201")[0]
    cdna = mtnd5_001.coding_sequence
    amino_acids, ends_with_stop_codon = translate_cdna(
        cdna,
        first_codon_is_start=True,
        translation_table=2)
    assert ends_with_stop_codon
    eq_(amino_acids, mtnd5_001.protein_sequence)


# Independent published TCAG-order strings, not derived from Isovar tables.
CODONS = ["".join(bases) for bases in product("TCAG", repeat=3)]
NCBI_TABLE_1 = "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"
NCBI_TABLE_2 = "FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIMMTTTTNNKKSS**VVVVAAAADDEEGGGG"


@pytest.mark.parametrize("table,assignments", [(1, NCBI_TABLE_1), (2, NCBI_TABLE_2)])
@pytest.mark.parametrize("index", range(64))
def test_every_nuclear_and_mitochondrial_codon(index, table, assignments):
    aa = assignments[index]
    expected = ("", True) if aa == "*" else (aa, False)
    assert translate_cdna(CODONS[index], translation_table=table) == expected
    if table == 1:
        # SV, CTA, transcript-frame and downstream Vaxrank consumers use this
        # public object instead of translate_cdna; they must decode identically.
        assert standard_genetic_code.translate(CODONS[index]) == expected
        assert translate_cdna(CODONS[index]) == expected


@pytest.mark.parametrize("stop", ["AGA", "AGG", "TAA", "TAG"])
def test_mitochondrial_termination_does_not_translate_downstream_codons(stop):
    assert translate_cdna("ATG" + stop + "GCT", translation_table=2) == ("M", True)


@pytest.mark.parametrize("start", ["ATT", "ATC", "ATA", "ATG", "GTG"])
def test_mitochondrial_initiation_does_not_change_elongation(start):
    internal = NCBI_TABLE_2[CODONS.index(start)]
    assert translate_cdna(start + start, translation_table=2, first_codon_is_start=True) == ("M" + internal, False)
    assert translate_cdna(start + start, translation_table=2, first_codon_is_start=False) == (internal * 2, False)


@pytest.mark.parametrize("start", ["ATG", "CTG", "TTG"])
def test_standard_initiators_match_the_ncbi_table(start):
    internal = NCBI_TABLE_1[CODONS.index(start)]
    assert standard_genetic_code.start_codons == {"ATG", "CTG", "TTG"}
    expected = ("M" + internal, False)
    assert standard_genetic_code.translate(start + start, first_codon_is_start=True) == expected
    assert translate_cdna(start + start, translation_table=1, first_codon_is_start=True) == expected
    assert translate_cdna(start + start, True) == expected


@pytest.mark.parametrize("start", ["GTG", "AGG", "ACG", "AAG", "ATC", "ATA", "ATT"])
def test_table_one_does_not_silently_expand_its_initiator_set(start):
    internal = NCBI_TABLE_1[CODONS.index(start)]
    assert translate_cdna(start + start, translation_table=1, first_codon_is_start=True) == (internal * 2, False)


def test_default_translation_uses_standard_nuclear_code():
    assert translate_cdna("ATATGAAGG") == ("I", True)
    assert translate_cdna("ATATGAAGG", translation_table=2) == ("MW", True)


def test_removed_mitochondrial_flag_is_rejected():
    with pytest.raises(TypeError, match="mitochondrial"):
        translate_cdna("ATGTGA", translation_table=2, mitochondrial=True)


def test_undeclared_annotated_code_keeps_standard_decoding():
    # Released ORFdata/Isovar exports and Vaxrank validation retain null source
    # metadata while using the producer's existing standard decoding convention.
    assert translate_cdna("ATGTTTTAA", translation_table=None) == ("MF", True)
    assert translate_cdna("ATATGAAGG", translation_table=None) == ("I", True)
    assert translate_cdna("GTGGTGTAA", translation_table=None, first_codon_is_start=True) == ("VV", True)


def test_declared_custom_start_does_not_mutate_the_shared_code():
    custom = standard_genetic_code.copy("source-declared-GTG", start_codons={"GTG"})
    assert custom.translate("GTGGTGTAA", first_codon_is_start=True) == ("MV", True)
    assert standard_genetic_code.translate("GTGGTGTAA", first_codon_is_start=True) == ("VV", True)
    assert translate_cdna("GTGGTGTAA", translation_table=1, first_codon_is_start=True) == ("VV", True)
