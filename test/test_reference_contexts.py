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
    reference_contexts_for_variants,
    variants_to_reference_contexts_dataframe,
    ReferenceContext,
)
from varcode import Variant, VariantCollection
from pyensembl import ensembl_grch38
from nose.tools import eq_

from testing_helpers import assert_equal_fields, load_vcf

def test_sequence_key_with_reading_frame_substitution_on_negative_strand():
    # replace second codon of TP53-001 with 'CCC'
    tp53_substitution = Variant(
        "17", 7676589, "CTC", "GGG", ensembl_grch38)
    variant_collection = VariantCollection([tp53_substitution])

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

    # first calling without a transcript ID white to see if we get back
    # multiple contexts
    reference_context_dict_many_transcripts = \
        reference_contexts_for_variants(
            variants=variant_collection,
            context_size=10,
            transcript_id_whitelist=None)

    assert len(reference_context_dict_many_transcripts) == 1, \
        "Dictionary should have only one variant but got %d keys" % (
            len(reference_context_dict_many_transcripts),)

    reference_contexts = reference_context_dict_many_transcripts[tp53_substitution]

    assert len(reference_contexts) > 1, \
        "Expected multiple reference contexts for %s but got %d: %s" % (
            tp53_substitution,
            len(reference_contexts),
            reference_contexts)

    reference_context_dict_single_transcript = \
        reference_contexts_for_variants(
            variants=variant_collection,
            context_size=10,
            transcript_id_whitelist={tp53_001.id})

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
    assert_equal_fields(result, expected)

def test_variants_to_reference_contexts_dataframe():
    variants = load_vcf("data/b16.f10/b16.vcf")
    assert len(variants) > 0
    df = variants_to_reference_contexts_dataframe(variants, context_size=10)
    print(df)
    groups = df.groupby(["chr", "pos", "ref", "alt"])
    # make sure we have at least one reference context for each
    # of the B16 coding variants
    eq_(len(groups), len(variants))


if __name__ == "__main__":
    test_sequence_key_with_reading_frame_substitution_on_negative_strand()
    test_variants_to_reference_contexts_dataframe()
