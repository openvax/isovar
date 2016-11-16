from nose.tools import eq_
from varcode import Variant
from isovar.variant_sequence_in_reading_frame import (
    compute_offset_to_first_complete_codon,
    match_variant_sequence_to_reference_context,
    VariantSequenceInReadingFrame,
)
from isovar.variant_sequences import VariantSequence
from isovar.reference_coding_sequence_key import ReferenceCodingSequenceKey
from isovar.reference_context import ReferenceContext
from isovar.allele_reads import AlleleRead


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

def make_inputs_for_tp53_201_variant(
        context_size=3,
        prefix=None):
    # TP53-201 is an isoform of TP53 which seems to lack untranslated
    # regions so the sequence is:
    # First exon: chr17 7,676,594 - 7,676,521
    # ATG|GAG|GAG|CCG|CAG|TCA|GAT...
    # -M-|-E-|-E-|-P-|-Q-|-S-|-D-

    # we're assuming a variant
    # chr17. 7,676,591 C>T which changes GAG (E) > AAG (K)
    variant = Variant("chr17", 7676591, "C", "T", "GRCh38")

    # TP53-201
    transcript = variant.ensembl.transcripts_by_name("TP53-201")[0]

    effect = variant.effect_on_transcript(transcript)

    eq_(effect.__class__.__name__, "Substitution")
    eq_(effect.aa_ref, "E")
    eq_(effect.aa_alt, "K")

    prefix = "ATG" if prefix is None else prefix
    alt = "A"
    suffix = "AGGAGCCGCAGTCAGAT"
    cdna_sequence = prefix + alt + suffix
    # variant sequence supported by two reads
    # one fully spanning the variant sequence
    # and another missing the first and last
    # nucleotides
    variant_sequence = VariantSequence(
        prefix=prefix,
        alt=alt,
        suffix=suffix,
        reads=[
            AlleleRead(
                prefix=prefix, allele=alt, suffix=suffix,
                name="full-overlap"),
            AlleleRead(
                prefix=prefix[1:], allele=alt, suffix=suffix[:-1],
                name="partial-overlap"),
        ])
    assert isinstance(variant_sequence, VariantSequence)

    reference_coding_sequence_key = ReferenceCodingSequenceKey.from_variant_and_transcript(
        variant=variant,
        transcript=transcript,
        context_size=3)
    assert isinstance(reference_coding_sequence_key, ReferenceCodingSequenceKey)

    reference_context = ReferenceContext.from_reference_coding_sequence_key(
        key=reference_coding_sequence_key,
        variant=variant,
        transcripts=[transcript])
    assert isinstance(reference_context, ReferenceContext)

    expected = VariantSequenceInReadingFrame(
        cdna_sequence=cdna_sequence,
        offset_to_first_complete_codon=0,
        variant_cdna_interval_start=3,
        variant_cdna_interval_end=4,
        reference_cdna_sequence_before_variant="ATG",
        number_mismatches=0)
    assert isinstance(expected, VariantSequenceInReadingFrame)

    return variant_sequence, reference_context, expected

def test_match_variant_sequence_to_reference_context_exact_match():
    variant_sequence, reference_context, expected = \
        make_inputs_for_tp53_201_variant(context_size=3)

    result = match_variant_sequence_to_reference_context(
        variant_sequence=variant_sequence,
        reference_context=reference_context,
        min_transcript_prefix_length=3,
        max_transcript_mismatches=0)
    eq_(expected, result)
