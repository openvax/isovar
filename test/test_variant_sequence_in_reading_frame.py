from nose.tools import eq_
from isovar.variant_sequence_in_reading_frame import (
    compute_offset_to_first_complete_codon,
    match_variant_sequence_to_reference_context,
    VariantSequenceInReadingFrame
)

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

def test_match_variant_sequence_to_reference_context_exact_match():
    # Coding sequence = --|ATG|CCC|TAG
    # first two nucleotides are untranslated region
    cdna_sequence = "GGATGCCCTAG"
    varseq_in_orf = VariantSequenceInReadingFrame(
        cdna_sequence=cdna_sequence,
        offset_to_first_complete_codon=2,
        variant_cdna_interval_start=5,
        variant_cdna_interval_end=6,
        reference_cdna_sequence_before_variant="GGATG",
        number_mismatches=0)

