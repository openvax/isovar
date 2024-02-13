from __future__ import print_function, division, absolute_import

from varcode import Variant
from nose.tools import eq_

from isovar.reference_sequence_key import ReferenceSequenceKey

from genomes_for_testing import grch38, grcm38


def test_sequence_key_for_variant_on_transcript_substitution():
    # rs769125639 is a simple T>A substitution in the 6th nucleotide of
    # BRCA2-001's 5' UTR
    brca2_variant_rs769125639 = Variant("13", 32315479, "T", "A", grch38)
    brca2_001 = grch38.transcripts_by_name("BRCA2-001")[0]
    # first 50 characters of BRCA2-001:
    #  "GGGCTTGTGGCGCGAGCTTCTGAAACTAGGCGGCAGAGGCGGAGCCGCTG"
    brca2_ref_seq = brca2_001.sequence[:50]
    eq_(brca2_ref_seq, "GGGCTTGTGGCGCGAGCTTCTGAAACTAGGCGGCAGAGGCGGAGCCGCTG")
    print(brca2_ref_seq)
    # get the 5 nucleotides before the variant and 10 nucleotides after
    sequence_key = ReferenceSequenceKey.from_variant_and_transcript(
        variant=brca2_variant_rs769125639,
        transcript=brca2_001,
        context_size=10)
    expected_sequence_key = ReferenceSequenceKey(
        strand="+",
        sequence_before_variant_locus=brca2_ref_seq[:5],
        sequence_at_variant_locus="T",
        sequence_after_variant_locus=brca2_ref_seq[6:16])
    eq_(sequence_key, expected_sequence_key)


def test_sequence_key_for_variant_on_transcript_deletion():
    # Delete the 6th nucleotide of BRCA2-001's 5' UTR
    brca2_variant_deletion = Variant("13", 32315479, "T", "", grch38)
    brca2_001 = grch38.transcripts_by_name("BRCA2-001")[0]
    # first 50 characters of BRCA2-001:
    #  "GGGCTTGTGGCGCGAGCTTCTGAAACTAGGCGGCAGAGGCGGAGCCGCTG"
    brca2_ref_seq = brca2_001.sequence[:50]
    eq_(brca2_ref_seq, "GGGCTTGTGGCGCGAGCTTCTGAAACTAGGCGGCAGAGGCGGAGCCGCTG")
    print(brca2_ref_seq)
    # get the 5 nucleotides before the variant and 10 nucleotides after
    sequence_key = ReferenceSequenceKey.from_variant_and_transcript(
        variant=brca2_variant_deletion,
        transcript=brca2_001,
        context_size=10)
    expected_sequence_key = ReferenceSequenceKey(
        strand="+",
        sequence_before_variant_locus=brca2_ref_seq[:5],
        sequence_at_variant_locus="T",
        sequence_after_variant_locus=brca2_ref_seq[6:16])
    eq_(sequence_key, expected_sequence_key)


def test_sequence_key_for_variant_on_transcript_insertion():
    # Insert 'CCC' after the 6th nucleotide of BRCA2-001's 5' UTR
    brca2_variant_insertion = Variant(
        "13", 32315479, "T", "TCCC", grch38)
    brca2_001 = grch38.transcripts_by_name("BRCA2-001")[0]
    # first 50 characters of BRCA2-001:
    #  "GGGCTTGTGGCGCGAGCTTCTGAAACTAGGCGGCAGAGGCGGAGCCGCTG"
    brca2_ref_seq = brca2_001.sequence[:50]
    eq_(brca2_ref_seq, "GGGCTTGTGGCGCGAGCTTCTGAAACTAGGCGGCAGAGGCGGAGCCGCTG")
    print(brca2_ref_seq)
    # get the 5 nucleotides before the variant and 10 nucleotides after
    sequence_key = ReferenceSequenceKey.from_variant_and_transcript(
        variant=brca2_variant_insertion,
        transcript=brca2_001,
        context_size=10)

    # expecting nothing at the variant locus since we're inserting between
    # two reference nucleotides
    expected_sequence_key = ReferenceSequenceKey(
        strand="+",
        sequence_before_variant_locus=brca2_ref_seq[:6],
        sequence_at_variant_locus="",
        sequence_after_variant_locus=brca2_ref_seq[6:16])
    eq_(sequence_key, expected_sequence_key)


def test_sequence_key_for_variant_on_transcript_substitution_reverse_strand():
    # Replace start codon of TP53-001 with 'CCC', however since this is on
    # reverse strand the variant becomes "CAT">"GGG"
    tp53_substitution = Variant("17", 7676592, "CAT", "GGG", grch38)
    tp53_001 = grch38.transcripts_by_name("TP53-001")[0]
    # Sequence of TP53 around start codon with 10 context nucleotides:
    # In [51]: t.sequence[190-10:190+13]
    # Out[51]: 'GGTCACTGCC_ATG_GAGGAGCCGC'
    eq_(tp53_001.sequence[190 - 10:190 + 13], "GGTCACTGCCATGGAGGAGCCGC")

    # get the 5 nucleotides before the variant and 10 nucleotides after
    sequence_key = ReferenceSequenceKey.from_variant_and_transcript(
        variant=tp53_substitution,
        transcript=tp53_001,
        context_size=10)

    expected_sequence_key = ReferenceSequenceKey(
        strand="-",
        sequence_before_variant_locus="GGTCACTGCC",
        sequence_at_variant_locus="ATG",
        sequence_after_variant_locus="GAGGAGCCGC")
    eq_(sequence_key, expected_sequence_key)


def test_sequence_key_for_variant_on_transcript_deletion_reverse_strand():
    # delete start codon of TP53-001, which in reverse complement means
    # deleting the sequence "CAT"
    tp53_deletion = Variant("17", 7676592, "CAT", "", grch38)
    tp53_001 = grch38.transcripts_by_name("TP53-001")[0]
    # Sequence of TP53 around start codon with 10 context nucleotides:
    # In [51]: t.sequence[190-10:190+13]
    # Out[51]: 'GGTCACTGCC_ATG_GAGGAGCCGC'
    eq_(tp53_001.sequence[190 - 10:190 + 13], "GGTCACTGCCATGGAGGAGCCGC")

    # get the 5 nucleotides before the variant and 10 nucleotides after
    sequence_key = ReferenceSequenceKey.from_variant_and_transcript(
        variant=tp53_deletion,
        transcript=tp53_001,
        context_size=10)

    expected_sequence_key = ReferenceSequenceKey(
        strand="-",
        sequence_before_variant_locus="GGTCACTGCC",
        sequence_at_variant_locus="ATG",
        sequence_after_variant_locus="GAGGAGCCGC")
    eq_(sequence_key, expected_sequence_key)


def test_sequence_key_for_variant_on_transcript_insertion_reverse_strand():
    # insert 'CCC' after start codon of TP53-001, which on the reverse
    # complement means inserting "GGG" between "CTC_CAT"
    tp53_insertion = Variant("17", 7676589, "CTC", "CTCGGG", grch38)
    tp53_001 = grch38.transcripts_by_name("TP53-001")[0]
    # Sequence of TP53 around start codon with 10 context nucleotides:
    # In [51]: t.sequence[190-10:190+13]
    # Out[51]: 'GGTCACTGCC_ATG_GAGGAGCCGC'
    eq_(tp53_001.sequence[190 - 10:190 + 13], "GGTCACTGCCATGGAGGAGCCGC")

    # The above gives us the cDNA sequence from the transcript, whereas the
    # reverse complement genomic sequence is:
    #    GCGGCTCCTC_CAT_GGCAGTGACC

    # get the 5 nucleotides before the variant and 10 nucleotides after
    sequence_key = ReferenceSequenceKey.from_variant_and_transcript(
        variant=tp53_insertion,
        transcript=tp53_001,
        context_size=10)

    expected_sequence_key = ReferenceSequenceKey(
        strand="-",
        sequence_before_variant_locus="CACTGCCATG",
        sequence_at_variant_locus="",
        sequence_after_variant_locus="GAGGAGCCGC")
    eq_(sequence_key, expected_sequence_key)


def test_reference_sequence_key_hash_and_equality_same_objects():
    rsk1 = ReferenceSequenceKey(
        strand="+",
        sequence_before_variant_locus="AAA",
        sequence_at_variant_locus="T",
        sequence_after_variant_locus="GGG")
    rsk2 = ReferenceSequenceKey(
        strand="+",
        sequence_before_variant_locus="AAA",
        sequence_at_variant_locus="T",
        sequence_after_variant_locus="GGG")
    eq_(rsk1, rsk2)
    eq_(hash(rsk1), hash(rsk2))
    eq_(str(rsk1), str(rsk2))
    eq_(repr(rsk1), repr(rsk2))


def test_reference_sequence_key_hash_and_equality_different_objects():
    rsk1 = ReferenceSequenceKey(
        strand="+",
        sequence_before_variant_locus="AAA",
        sequence_at_variant_locus="T",
        sequence_after_variant_locus="GGG")
    rsk_different_strand = ReferenceSequenceKey(
        strand="-",
        sequence_before_variant_locus="AAA",
        sequence_at_variant_locus="T",
        sequence_after_variant_locus="GGG")
    assert rsk1 != rsk_different_strand
    assert hash(rsk1) != hash(rsk_different_strand)
    assert str(rsk1) != str(rsk_different_strand)
    assert repr(rsk1) != repr(rsk_different_strand)


def test_reference_sequence_key_from_weird_deletion():
    # variant reads into the intron; want to make sure isovar skips over such cases

    variant = Variant("11", 106262686, "GTGAAGG", "", grcm38)
    transcript = grcm38.transcript_by_id("ENSMUST00000021049")
    sequence_key = ReferenceSequenceKey.from_variant_and_transcript(
        variant=variant,
        transcript=transcript,
        context_size=10)
    assert sequence_key is None, '%s\n%s' % (sequence_key, transcript)
