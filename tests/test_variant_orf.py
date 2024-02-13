from __future__ import print_function, division, absolute_import

from nose.tools import eq_
from varcode import Variant

from isovar.variant_orf import (
    compute_offset_to_first_complete_codon,
    VariantORF,
)
from isovar.variant_orf_helpers import match_variant_sequence_to_reference_context
from isovar.variant_sequence import VariantSequence
from isovar.reference_coding_sequence_key import ReferenceCodingSequenceKey
from isovar.reference_context import ReferenceContext
from isovar.allele_read import AlleleRead
from isovar.dna import reverse_complement_dna

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
        cdna_prefix="ATG",
        cdna_suffix="AGGAGCCGCAGTCAGAT",
        n_bad_nucleotides_at_start=0,
        mismatches_before_variant=0,
        mismatches_after_variant=14,  # the read is that much longer than the reference (17 vs 3)
        reference_context_size=3):
    """
    Parameters
    ----------
    cdna_prefix : str
        Transcript nucleotides before the variant that we're pretending
        got detected from RNA-seq reads.

    cdna_suffix : str
        Transcript nucleotides after the variant that we're pretending
        got detected from RNA-seq reads.

    n_bad_nucleotides_at_start : int
        Number of nucleotides we expect to get trimmed from the
        beginning of the variant sequence while matching to a reference context.

    mismatches_before_variant : int
        Expected number of nucleotide mismatches in the result before
        the variant locus.

    reference_context_size : int
        Number of nucleotides before the variant locus to try matching
        against a reference transcript.
    """
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

    cdna_alt = "A"

    # genomic DNA is the reverse complement of the cDNA
    # for TP53-001 since it's on the negative strand
    gdna_prefix = reverse_complement_dna(cdna_suffix)
    gdna_alt = reverse_complement_dna(cdna_alt)
    gdna_suffix = reverse_complement_dna(cdna_prefix)

    # variant sequence supported by two reads
    # one fully spanning the variant sequence
    # and another missing the last nucleotide
    fully_overlapping_read = AlleleRead(
        prefix=gdna_prefix,
        allele=gdna_alt,
        suffix=gdna_suffix,
        name="full-overlap")
    # testing the prefix and allele to make sure they have the expected
    # TP53-201 sequence but the suffix might change depending on what's
    # passed in as cdna_prefix
    if cdna_suffix == "AGGAGCCGCAGTCAGAT":
        eq_(fully_overlapping_read.prefix, "ATCTGACTGCGGCTCCT")
    eq_(fully_overlapping_read.allele, "T")

    partially_overlapping_read = AlleleRead(
        prefix=gdna_prefix,
        allele=gdna_alt,
        suffix=gdna_suffix[:-1],
        name="partial-overlap")
    if cdna_suffix == "AGGAGCCGCAGTCAGAT":
        eq_(partially_overlapping_read.prefix, "ATCTGACTGCGGCTCCT")
    eq_(partially_overlapping_read.allele, "T")

    variant_sequence = VariantSequence(
        prefix=gdna_prefix,
        alt=gdna_alt,
        suffix=gdna_suffix,
        reads=[fully_overlapping_read, partially_overlapping_read])
    assert isinstance(variant_sequence, VariantSequence)

    prefix_length = len(cdna_prefix) - n_bad_nucleotides_at_start

    reference_coding_sequence_key = ReferenceCodingSequenceKey.from_variant_and_transcript(
        variant=variant,
        transcript=transcript,
        context_size=reference_context_size)
    assert isinstance(reference_coding_sequence_key, ReferenceCodingSequenceKey)

    reference_context = ReferenceContext.from_reference_coding_sequence_key(
        key=reference_coding_sequence_key,
        variant=variant,
        transcripts=[transcript])
    assert isinstance(reference_context, ReferenceContext)

    expected = VariantORF(
        cdna_sequence=cdna_prefix[-prefix_length:] + cdna_alt + cdna_suffix,
        offset_to_first_complete_codon=prefix_length % 3,
        variant_cdna_interval_start=prefix_length,
        variant_cdna_interval_end=prefix_length + 1,
        reference_cdna_sequence_before_variant="ATG"[-prefix_length:],
        reference_cdna_sequence_after_variant="AGGAGCCGCAGTCAGAT"[:reference_context_size],
        num_mismatches_before_variant=mismatches_before_variant,
        num_mismatches_after_variant=mismatches_after_variant)
    assert isinstance(expected, VariantORF)

    return variant_sequence, reference_context, expected


def test_match_variant_sequence_to_reference_context_exact_match():
    # Variant sequence is exact match for beginning of TP53-201 transcript
    variant_sequence, reference_context, expected = \
        make_inputs_for_tp53_201_variant()

    result = match_variant_sequence_to_reference_context(
        variant_sequence=variant_sequence,
        reference_context=reference_context,
        min_transcript_prefix_length=3,
        max_transcript_mismatches=0)
    eq_(expected, result)


def test_match_variant_sequence_to_reference_context_not_enough_prefix():
    # Variant sequence missing first nucleotide of start codon
    # ("TG" instead of "ATG") and the variant occurrs immediately after
    # the start codon. Since the min_transcript_prefix_length is 3 in
    # this case we expect the match function to return None
    variant_sequence, reference_context, _ = \
        make_inputs_for_tp53_201_variant(
            cdna_prefix="TG",
            reference_context_size=2)

    result = match_variant_sequence_to_reference_context(
        variant_sequence=variant_sequence,
        reference_context=reference_context,
        min_transcript_prefix_length=3,
        max_transcript_mismatches=0)
    eq_(result, None)


def test_match_variant_sequence_to_reference_context_trim_1_bad_nucleotide():
    # Variant sequence has an extra nucleotide at the beginning which is
    # supported by only 1 read, whereas the correct sequence is supported by
    # 2 reads. If we allow > 1 "attempt" in the match function then it will
    # trim off the extra "G" and correctly match against the TP53-201
    # transcript sequence.

    variant_sequence, reference_context, expected = \
        make_inputs_for_tp53_201_variant(
            cdna_prefix="GATG",
            n_bad_nucleotides_at_start=1)

    result = match_variant_sequence_to_reference_context(
        variant_sequence=variant_sequence,
        reference_context=reference_context,
        min_transcript_prefix_length=3,
        max_transcript_mismatches=0,
        max_trimming_attempts=1)
    eq_(expected, result)


def test_match_variant_sequence_to_reference_context_ignore_extra_prefix():
    # There are three "extra" nucleotides at the start but since we are
    # only using reference context size of 3 then this sequence will
    # match.
    variant_sequence, reference_context, expected = \
        make_inputs_for_tp53_201_variant(
            cdna_prefix="GGGATG",
            n_bad_nucleotides_at_start=3,
            reference_context_size=3)

    result = match_variant_sequence_to_reference_context(
        variant_sequence=variant_sequence,
        reference_context=reference_context,
        min_transcript_prefix_length=3,
        max_transcript_mismatches=0,
        max_trimming_attempts=0)
    eq_(expected, result)
    # make sure that the "GGG" codon got ignored since translation
    # should start at the "ATG" after it
    eq_(result.cdna_sequence[:3], "ATG")


def test_match_variant_sequence_to_reference_context_bad_start_nucleotide_no_trimming():
    # matching should fail if no mismatches are allowed and no trimming rounds
    # are allowed
    variant_sequence, reference_context, _ = \
        make_inputs_for_tp53_201_variant(
            cdna_prefix="CTG",
            n_bad_nucleotides_at_start=1)

    result = match_variant_sequence_to_reference_context(
        variant_sequence=variant_sequence,
        reference_context=reference_context,
        min_transcript_prefix_length=2,
        max_transcript_mismatches=0,
        max_trimming_attempts=0)
    eq_(None, result)


def test_match_variant_sequence_to_reference_context_bad_start_nucleotide_trimming():
    # match should succeed if 1 round of trimming is allowed
    variant_sequence, reference_context, expected = \
        make_inputs_for_tp53_201_variant(
            cdna_prefix="CTG",
            n_bad_nucleotides_at_start=1)
    result = match_variant_sequence_to_reference_context(
        variant_sequence=variant_sequence,
        reference_context=reference_context,
        min_transcript_prefix_length=2,
        max_transcript_mismatches=0,
        max_trimming_attempts=1)
    eq_(expected, result)


def test_match_variant_sequence_to_reference_context_bad_start_nucleotide_allow_mismatch():
    # match should succeed if 1 mismatch is allowed
    variant_sequence, reference_context, expected = \
        make_inputs_for_tp53_201_variant(
            cdna_prefix="CTG",
            n_bad_nucleotides_at_start=0,
            mismatches_before_variant=1)
    result = match_variant_sequence_to_reference_context(
        variant_sequence=variant_sequence,
        reference_context=reference_context,
        min_transcript_prefix_length=3,
        max_transcript_mismatches=1,
        max_trimming_attempts=0)
    eq_(expected, result)


def test_match_variant_sequence_to_reference_context_include_mismatches_after_variant():
    variant_sequence, reference_context, expected = \
        make_inputs_for_tp53_201_variant(
            cdna_suffix="AGAAGCCGCAGTCAGAT",  # too long and also one mismatch: G>A in 3rd char
            mismatches_after_variant=15)

    result = match_variant_sequence_to_reference_context(
        variant_sequence=variant_sequence,
        reference_context=reference_context,
        min_transcript_prefix_length=3,
        max_transcript_mismatches=0,
        count_mismatches_after_variant=False)
    # should have a result, since we're not counting mismatches after the variant
    eq_(expected, result)

    # now say we want to count mismatches after the variant - expect no result
    result = match_variant_sequence_to_reference_context(
        variant_sequence=variant_sequence,
        reference_context=reference_context,
        min_transcript_prefix_length=3,
        max_transcript_mismatches=0,
        count_mismatches_after_variant=True)
    eq_(None, result)
