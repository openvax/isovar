# Copyright (c) 2016. Mount Sinai School of Medicine
#
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

from __future__ import print_function, division, absolute_import

from isovar.reference_context import (
    reading_frame_to_offset,
    sequence_key_with_reading_frame_for_variant_on_transcript,
    SequenceKeyWithReadingFrame,
)
from varcode import Variant
from pyensembl import ensembl_grch38
from nose.tools import eq_

from testing_helpers import assert_equal_fields

def test_reading_frame_to_offset():
    eq_(reading_frame_to_offset(0), 0)
    eq_(reading_frame_to_offset(1), 2)
    eq_(reading_frame_to_offset(2), 1)


def test_sequence_key_with_reading_frame_substitution_with_five_prime_utr():
    # Replace second codon of TP53-001 with 'CCC', the surrounding context
    # includes nucleotides from the 5' UTR. Since TP53 is on the negative
    # strand we have to take the reverse complement of the variant which turns
    # it into CTC>GGG
    tp53_substitution = Variant(
        "17", 7676589, "CTC", "GGG", ensembl_grch38)
    tp53_001 = ensembl_grch38.transcripts_by_name("TP53-001")[0]

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
    result = \
        sequence_key_with_reading_frame_for_variant_on_transcript(
            variant=tp53_substitution,
            transcript=tp53_001,
            context_size=10)
    expected = SequenceKeyWithReadingFrame(
        strand="-",
        sequence_before_variant_locus="CACTGCCATG",
        sequence_at_variant_locus="GAG",
        sequence_after_variant_locus="GAGCCGCAGT",
        offset_to_first_complete_codon=7,
        contains_start_codon=True,
        overlaps_start_codon=True,
        contains_five_prime_utr=True,
        amino_acids_before_variant="M")
    assert_equal_fields(result, expected)

def test_sequence_key_with_reading_frame_deletion_with_five_prime_utr():
    # Delete second codon of TP53-001, the surrounding context
    # includes nucleotides from the 5' UTR. Since TP53 is on the negative
    # strand we have to take the reverse complement of the variant which turns
    # it into 'CTC'>''
    tp53_deletion = Variant(
        "17", 7676589, "CTC", "", ensembl_grch38)
    tp53_001 = ensembl_grch38.transcripts_by_name("TP53-001")[0]

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

    result = \
        sequence_key_with_reading_frame_for_variant_on_transcript(
            variant=tp53_deletion,
            transcript=tp53_001,
            context_size=10)
    expected = SequenceKeyWithReadingFrame(
        strand="-",
        sequence_before_variant_locus="CACTGCCATG",
        sequence_at_variant_locus="GAG",
        sequence_after_variant_locus="GAGCCGCAGT",
        offset_to_first_complete_codon=7,
        contains_start_codon=True,
        overlaps_start_codon=True,
        contains_five_prime_utr=True,
        amino_acids_before_variant="M")
    assert_equal_fields(result, expected)


def test_sequence_key_with_reading_frame_insertion():
    # Insert nucleotide "T" after second codon of TP53-001, the
    # surrounding context includes nucleotides from the 5' UTR. Since TP53 is on
    # the negative strand we have to take the reverse complement of the variant
    # which turns it into 'CTC'>'CTCA'
    tp53_insertion = Variant(
        "17", 7676586, "CTC", "CTCA", ensembl_grch38)

    tp53_001 = ensembl_grch38.transcripts_by_name("TP53-001")[0]
    # Sequence of TP53 around boundary of 2nd/3rd codons
    # with 10 context nucleotides:
    #   last 4 nt of 5' UTR: TGCC
    #   start codon: ATG (translates to M)
    #   2nd codon: GAG (translates to E)
    #   <---- insertion variant occurs between these two codons
    #   3rd codon: GAG
    #   4th codon: CCG
    #   5th codon:  CAG
    #   first nt of 6th codon: T

    result = \
        sequence_key_with_reading_frame_for_variant_on_transcript(
            variant=tp53_insertion,
            transcript=tp53_001,
            context_size=10)

    expected = SequenceKeyWithReadingFrame(
        strand="-",
        sequence_before_variant_locus="TGCCATGGAG",
        sequence_at_variant_locus="",
        sequence_after_variant_locus="GAGCCGCAGT",
        offset_to_first_complete_codon=4,
        contains_start_codon=True,
        overlaps_start_codon=True,
        contains_five_prime_utr=True,
        amino_acids_before_variant="ME")
    assert_equal_fields(result, expected)

def test_sequence_key_with_reading_frame_insertion_inside_start_codon():
    # insert nucleotide "C" in the middle of the start codon of TP53-001,
    # keeping only 1 nucleotide of context. In the reverse complement this
    # becomes 'T'>'TG'
    tp53_insertion = Variant(
        "17", 7676592, "T", "TG", ensembl_grch38)

    tp53_001 = ensembl_grch38.transcripts_by_name("TP53-001")[0]

    result = \
        sequence_key_with_reading_frame_for_variant_on_transcript(
            variant=tp53_insertion,
            transcript=tp53_001,
            context_size=1)
    assert result is None, "Expected result to be None when variant affects start codon"

def test_sequence_key_with_reading_frame_insertion_before_start_codon():
    # insert nucleotide "T" before of the start codon of TP53-001,
    tp53_insertion = Variant("17", 7676593, "C", "CT", ensembl_grch38)

    tp53_001 = ensembl_grch38.transcripts_by_name("TP53-001")[0]

    result = \
        sequence_key_with_reading_frame_for_variant_on_transcript(
            variant=tp53_insertion,
            transcript=tp53_001,
            context_size=1)
    assert result is None, "Expected result to be None when variant before start codon"


def test_sequence_key_with_reading_frame_insertion_context_6nt_contains_start():
    # Insert nucleotide "T" after second codon of TP53-001,
    # but in this test we're going to only keep enough context to see
    # the start codon but none of the 5' UTR. In the reverse complement this
    # variant becomes CTC>CTCA
    tp53_insertion = Variant(
        "17", 7676586, "CTC", "CTCA", ensembl_grch38)

    tp53_001 = ensembl_grch38.transcripts_by_name("TP53-001")[0]
    # Sequence of TP53 around boundary of 2nd/3rd codons
    # with 6 context nucleotides:
    #   start codon: ATG (translates to M)
    #   2nd codon: GAG (translates to E)
    #   <---- insertion variant occurs between these two codons
    #   3rd codon: GAG
    #   4th codon: CCG

    result = \
        sequence_key_with_reading_frame_for_variant_on_transcript(
            variant=tp53_insertion,
            transcript=tp53_001,
            context_size=6)

    expected = SequenceKeyWithReadingFrame(
        strand="-",
        sequence_before_variant_locus="ATGGAG",
        sequence_at_variant_locus="",
        sequence_after_variant_locus="GAGCCG",
        offset_to_first_complete_codon=0,
        contains_start_codon=True,
        overlaps_start_codon=True,
        contains_five_prime_utr=False,
        amino_acids_before_variant="ME")
    assert_equal_fields(result, expected)


def test_sequence_key_with_reading_frame_insertion_context_5nt_overlaps_start():
    # Insert nucleotide "T" after second codon of TP53-001,
    # but in this test we're going to only keep enough context to see
    # a part of the start codon, thus the result shouldn't "contain"
    # the start codon but does "overlap" it. In the reverse complement
    # this variant becomes CTC>CTCA
    tp53_insertion = Variant(
        "17", 7676586, "CTC", "CTCA", ensembl_grch38)

    tp53_001 = ensembl_grch38.transcripts_by_name("TP53-001")[0]
    # Sequence of TP53 around boundary of 2nd/3rd codons
    # with 6 context nucleotides:
    #   last two nt of start codon: TG
    #   2nd codon: GAG (translates to E)
    #   <---- insertion variant occurs between these two codons
    #   3rd codon: GAG
    #   first two nt of 4th codon: CC

    result = \
        sequence_key_with_reading_frame_for_variant_on_transcript(
            variant=tp53_insertion,
            transcript=tp53_001,
            context_size=5)

    expected = SequenceKeyWithReadingFrame(
        strand="-",
        sequence_before_variant_locus="TGGAG",
        sequence_at_variant_locus="",
        sequence_after_variant_locus="GAGCC",
        offset_to_first_complete_codon=2,
        contains_start_codon=False,
        overlaps_start_codon=True,
        contains_five_prime_utr=False,
        amino_acids_before_variant="E")
    assert_equal_fields(result, expected)


def test_sequence_key_with_reading_frame_insertion_context_3nt_no_start():
    # Insert nucleotide "T" after second codon of TP53-001,
    # but in this test we're going to only keep enough context to see
    # the second codon (and no nucleotides from the start). In the reverse
    # complement this variant becomes CTC>CTCA.

    tp53_insertion = Variant(
        "17", 7676586, "CTC", "CTCA", ensembl_grch38)

    tp53_001 = ensembl_grch38.transcripts_by_name("TP53-001")[0]
    # Sequence of TP53 around boundary of 2nd/3rd codons
    # with 6 context nucleotides:
    #   2nd codon: GAG (translates to E)
    #   <---- insertion variant occurs between these two codons
    #   3rd codon: GAG

    result = \
        sequence_key_with_reading_frame_for_variant_on_transcript(
            variant=tp53_insertion,
            transcript=tp53_001,
            context_size=3)

    expected = SequenceKeyWithReadingFrame(
        strand="-",
        sequence_before_variant_locus="GAG",
        sequence_at_variant_locus="",
        sequence_after_variant_locus="GAG",
        offset_to_first_complete_codon=0,
        contains_start_codon=False,
        overlaps_start_codon=False,
        contains_five_prime_utr=False,
        amino_acids_before_variant="E")
    assert_equal_fields(result, expected)
