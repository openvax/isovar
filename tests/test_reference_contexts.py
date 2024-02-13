from __future__ import print_function, division, absolute_import

from varcode import Variant, VariantCollection

from nose.tools import eq_

from isovar.reference_context import ReferenceContext
from isovar.reference_context_helpers import reference_contexts_generator

from isovar.dataframe_helpers import variants_to_reference_contexts_dataframe

from testing_helpers import load_vcf
from genomes_for_testing import grch38


def test_sequence_key_with_reading_frame_substitution_on_negative_strand():
    # replace second codon of TP53-001 with 'CCC'
    tp53_substitution = Variant(
        "17", 7676589, "CTC", "GGG", grch38)
    variant_collection = VariantCollection([tp53_substitution])

    tp53_001 = grch38.transcripts_by_name("TP53-001")[0]

    # Sequence of TP53 around second codon with 10 context nucleotides:
    # In [51]: t.sequence[193-10:193+13]
    # Out[51]: 'CACTGCCATGGAGGAGCCGCAGT'
    # Which can be split into the following parts:
    #  last 7 nt of 5' UTR: CACTGCC
    #  start codon: ATG (translates to M)
    #  2nd codon: GAG    <---- variant occurs here
    #  3rd codon: GAG
    #  4th codon: CCG
    #  5th codon:  CAG
    #  first nt of 6th codon: T

    # first calling without a transcript ID white to see if we get back
    # multiple contexts
    reference_contexts_gen = \
        reference_contexts_generator(
            variants=variant_collection,
            context_size=10,
            transcript_id_whitelist=None)

    reference_contexts_dict = dict(reference_contexts_gen)

    assert len(reference_contexts_dict) == 1, \
        "Dictionary should have only one variant but got %d keys" % (
            len(reference_contexts_dict),)

    reference_contexts = reference_contexts_dict[tp53_substitution]

    assert len(reference_contexts) > 1, \
        "Expected multiple reference contexts for %s but got %d: %s" % (
            tp53_substitution,
            len(reference_contexts),
            reference_contexts)

    reference_context_dict_single_transcript = \
        dict(reference_contexts_generator(
            variants=variant_collection,
            context_size=10,
            transcript_id_whitelist={tp53_001.id}))

    # still only expect one variant key
    eq_(len(reference_context_dict_single_transcript), 1)

    result_list = reference_context_dict_single_transcript[tp53_substitution]

    # since we limited the transcript ID whitelist, we only expect a single
    # reference context in the result
    eq_(len(result_list), 1)

    result = result_list[0]

    expected = ReferenceContext(
        strand="-",
        sequence_before_variant_locus="CACTGCCATG",
        sequence_at_variant_locus="GAG",
        sequence_after_variant_locus="GAGCCGCAGT",
        offset_to_first_complete_codon=7,
        contains_start_codon=True,
        overlaps_start_codon=True,
        contains_five_prime_utr=True,
        amino_acids_before_variant="M",
        variant=tp53_substitution,
        transcripts=[tp53_001])
    eq_(result, expected)


def test_variants_to_reference_contexts_dataframe():
    variants = load_vcf("data/b16.f10/b16.vcf")
    assert len(variants) > 0
    gen = reference_contexts_generator(variants, context_size=10)
    df = variants_to_reference_contexts_dataframe(gen)
    print(df)
    groups = df.groupby(["chr", "pos", "ref", "alt"])
    # make sure we have at least one reference context for each
    # of the B16 coding variants
    eq_(len(groups), len(variants))


def test_reference_context_hash_and_equality_same_object():
    rc1 = ReferenceContext(
        strand="-",
        sequence_before_variant_locus="C" * 7 + "ATG",
        sequence_at_variant_locus="T",
        sequence_after_variant_locus="TT",
        offset_to_first_complete_codon=7,
        contains_start_codon=True,
        overlaps_start_codon=True,
        contains_five_prime_utr=True,
        amino_acids_before_variant="M",
        variant=None,
        transcripts=[])
    rc2 = ReferenceContext(
        strand="-",
        sequence_before_variant_locus="C" * 7 + "ATG",
        sequence_at_variant_locus="T",
        sequence_after_variant_locus="TT",
        offset_to_first_complete_codon=7,
        contains_start_codon=True,
        overlaps_start_codon=True,
        contains_five_prime_utr=True,
        amino_acids_before_variant="M",
        variant=None,
        transcripts=[])

    eq_(rc1, rc2)
    eq_(str(rc1), str(rc2))
    eq_(repr(rc1), repr(rc2))
    eq_(hash(rc1), hash(rc2))


def test_reference_context_hash_and_equality_different_objects():
    rc1 = ReferenceContext(
        strand="-",
        sequence_before_variant_locus="C" * 7 + "ATG",
        sequence_at_variant_locus="T",
        sequence_after_variant_locus="TT",
        offset_to_first_complete_codon=7,
        contains_start_codon=True,
        overlaps_start_codon=True,
        contains_five_prime_utr=True,
        amino_acids_before_variant="M",
        variant=None,
        transcripts=[])

    rc_different_strand = ReferenceContext(
        strand="+",
        sequence_before_variant_locus="C" * 7 + "ATG",
        sequence_at_variant_locus="T",
        sequence_after_variant_locus="TT",
        offset_to_first_complete_codon=7,
        contains_start_codon=True,
        overlaps_start_codon=True,
        contains_five_prime_utr=True,
        amino_acids_before_variant="M",
        variant=None,
        transcripts=[])

    assert rc1 != rc_different_strand, \
        "Expected %s != %s" % (rc1, rc_different_strand)
    assert str(rc1) != str(rc_different_strand), \
        "Expected __str__ '%s' != '%s'" % (
            rc1, rc_different_strand)
    assert repr(rc1) != repr(rc_different_strand), \
        "Expected __repr__ '%s' != '%s'" % (
            rc1, rc_different_strand)

    assert hash(rc1) != hash(rc_different_strand), \
        "Expected hash(%s) != hash(%s) (%d vs. %d)" % (
            rc1,
            rc_different_strand,
            hash(rc1),
            hash(rc_different_strand))
