from isovar.translation_helpers import (
    compute_offset_to_first_complete_codon,
    translate_cdna)
from nose.tools import eq_
from pyensembl import ensembl_grch38

def test_compute_offset_to_first_complete_codon_no_trimming():
    # if nothing gets trimmed from the reference sequence, then
    # the offset to the first codon shouldn't change
    eq_(
        compute_offset_to_first_complete_codon(
            offset_to_first_complete_reference_codon=0,
            n_trimmed_from_reference_sequence=0),
        0)
    eq_(
        compute_offset_to_first_complete_codon(
            offset_to_first_complete_reference_codon=5,
            n_trimmed_from_reference_sequence=0),
        5)

def test_compute_offset_to_first_complete_codon_trimming_before_codon():
    # if the number of reference bases trimmed from the reference sequence
    # occurs before the reference codon, then it should decrease the
    # offset by the amount trimmed
    eq_(
        compute_offset_to_first_complete_codon(
            offset_to_first_complete_reference_codon=7,
            n_trimmed_from_reference_sequence=2),
        5)
    eq_(
        compute_offset_to_first_complete_codon(
            offset_to_first_complete_reference_codon=7,
            n_trimmed_from_reference_sequence=7),
        0)

def test_compute_offset_to_first_complete_codon_trimming_after_codon():
    # if the number of reference bases trimmed from the reference sequence
    # occurs after the reference codon, then it needs to be rounded up the
    # next multiple of three
    eq_(
        compute_offset_to_first_complete_codon(
            offset_to_first_complete_reference_codon=7,
            n_trimmed_from_reference_sequence=8),
        2)
    eq_(
        compute_offset_to_first_complete_codon(
            offset_to_first_complete_reference_codon=7,
            n_trimmed_from_reference_sequence=10),
        0)

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
    amino_acids = translate_cdna(cdna, first_codon_is_start=True)
    eq_(amino_acids, tp53_001.protein_sequence)
