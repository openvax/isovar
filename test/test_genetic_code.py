from isovar.genetic_code import translate_cdna
from nose.tools import eq_
from pyensembl import ensembl_grch38


def test_translate_cdna_no_stop_codon():
    eq_(translate_cdna("ATGATG", first_codon_is_start=False), ("MM", False))

def test_translate_cdna_stop_codon():
    eq_(translate_cdna("ATGATGTAG", first_codon_is_start=False), ("MM", True))

def test_translate_cdna_alternate_CTG_start():
    eq_(translate_cdna("CTGCTG", first_codon_is_start=True), ("ML", False))

def test_translate_cdna_CTG_after_start():
    eq_(translate_cdna("CTGCTG", first_codon_is_start=False), ("LL", False))

def test_TP53_translation_from_cdna():
    tp53_001 = ensembl_grch38.transcripts_by_name("TP53-001")[0]
    cdna = tp53_001.coding_sequence
    amino_acids, ends_with_stop_codon = translate_cdna(cdna, first_codon_is_start=True)
    assert ends_with_stop_codon
    eq_(amino_acids, tp53_001.protein_sequence)

def test_mitochondrial_MTND5_translation_from_cdna():
    mtnd5_001 = ensembl_grch38.transcripts_by_name("MT-ND5-201")[0]
    cdna = mtnd5_001.coding_sequence
    amino_acids, ends_with_stop_codon = translate_cdna(
        cdna,
        first_codon_is_start=True,
        mitochondrial=True)
    assert ends_with_stop_codon
    eq_(amino_acids, mtnd5_001.protein_sequence)
