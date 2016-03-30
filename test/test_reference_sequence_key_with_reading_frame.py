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

def test_reading_frame_to_offset():
    eq_(reading_frame_to_offset(0), 0)
    eq_(reading_frame_to_offset(1), 2)
    eq_(reading_frame_to_offset(2), 1)

def _check_equal_fields(result, expected):
    """
    Assert that fields of two SequenceKeyWithReadingFrame objects have
    same field values.
    """
    for field in SequenceKeyWithReadingFrame._fields:
        result_value = getattr(result, field)
        expected_value = getattr(expected, field)
        assert result_value == expected_value, \
            "Wrong value for '%s', expected %s but got %s" % (
                field,
                expected_value,
                result_value)


def test_sequence_key_with_reading_frame_substitution_with_five_prime_utr():
    # replace second codon of TP53-001 with 'CCC'
    tp53_substitution = Variant(
        "17", 7676589, "GAG", "CCC", ensembl_grch38)
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
    _check_equal_fields(result, expected)

def test_sequence_key_with_reading_frame_deletion_with_five_prime_utr():
    # delete second codon of TP53-001
    tp53_deletion = Variant(
        "17", 7676589, "GAG", "", ensembl_grch38)
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
    _check_equal_fields(result, expected)


def test_sequence_key_with_reading_frame_insertion_with_five_prime_utr():
    # insert nucleotide "T" after second codon of TP53-001
    tp53_insertion = Variant(
        "17", 7676586, "GAG", "GAGT", ensembl_grch38)

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
    _check_equal_fields(result, expected)
