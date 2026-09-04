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

import pytest
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
from isovar.protein_sequence_creator import ProteinSequenceCreator
from isovar.variant_sequence_creator import VariantSequenceCreator
from isovar.variant_sequence_helpers import filter_variant_sequences

from .common import eq_ 
from .genomes_for_testing import grch38

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


def test_filter_variant_sequences_defers_reference_compatibility():
    """
    A longer assembly can be unusable even when it has greater read support.

    Reference compatibility is unavailable during variant-sequence filtering,
    so the shorter compatible candidate must reach the translation stage.
    """
    base_sequence, reference_context, _ = make_inputs_for_tp53_201_variant()

    def with_support(suffix, n_fragments, name_prefix):
        reads = [
            AlleleRead(
                prefix=base_sequence.prefix,
                allele=base_sequence.alt,
                suffix=suffix,
                name="%s_%d" % (name_prefix, i))
            for i in range(n_fragments)
        ]
        return VariantSequence(
            prefix=base_sequence.prefix,
            alt=base_sequence.alt,
            suffix=suffix,
            reads=reads)

    compatible_short = with_support(
        suffix=base_sequence.suffix,
        n_fragments=3,
        name_prefix="compatible")
    # TP53-201 is on the negative strand. Prepending GGG to the genomic
    # suffix makes the three cDNA bases before the variant CCC instead of ATG,
    # exceeding max_transcript_mismatches=2 after length normalization.
    incompatible_long = with_support(
        suffix="GGG" + base_sequence.suffix,
        n_fragments=4,
        name_prefix="incompatible")

    filtered = filter_variant_sequences(
        variant_sequences=[compatible_short, incompatible_long],
        preferred_sequence_length=len(incompatible_long),
        min_variant_sequence_coverage=2)

    creator = ProteinSequenceCreator(
        protein_sequence_length=20,
        min_transcript_prefix_length=3,
        max_transcript_mismatches=2)
    translations = creator.all_pairs_translations(
        variant_sequences=filtered,
        reference_contexts=[reference_context])

    eq_(set(filtered), {compatible_short, incompatible_long})
    eq_(len(translations), 1)
    eq_(translations[0].untrimmed_variant_sequence, compatible_short)


@pytest.mark.parametrize("strand", ["+", "-"])
@pytest.mark.parametrize("variant_sequence_assembly", [True, False])
def test_nested_assemblies_reach_reference_matching_independently(
        strand,
        variant_sequence_assembly):
    """Regression test for #198 through sequence creation and translation."""
    variant = Variant("1", 100, "G", "C", "GRCh38")
    cdna_prefixes = ("AAA", "GGGAAA", "CCCGGGAAA", "TTTCCCGGGAAA")
    if strand == "+":
        genomic_prefixes = cdna_prefixes
        genomic_suffixes = ("A" * 8,) * len(cdna_prefixes)
    else:
        genomic_prefixes = ("T" * 8,) * len(cdna_prefixes)
        genomic_suffixes = tuple(
            reverse_complement_dna(prefix) for prefix in cdna_prefixes)
    reads = [
        AlleleRead(
            prefix=prefix,
            allele=variant.alt,
            suffix=suffix,
            name="%s_%s_%d" % (strand, cdna_prefix, i))
        for cdna_prefix, prefix, suffix in zip(
            cdna_prefixes, genomic_prefixes, genomic_suffixes)
        for i in range(2)
    ]
    reference_context = ReferenceContext(
        strand=strand,
        sequence_before_variant_locus="A" * len(cdna_prefixes[-1]),
        sequence_at_variant_locus=(
            variant.ref if strand == "+"
            else reverse_complement_dna(variant.ref)),
        sequence_after_variant_locus="A" * 8,
        offset_to_first_complete_codon=0,
        contains_start_codon=False,
        overlaps_start_codon=False,
        contains_five_prime_utr=False,
        amino_acids_before_variant="",
        variant=variant,
        transcripts=())

    sequence_creator = VariantSequenceCreator(
        min_variant_sequence_coverage=2,
        # VariantSequenceCreator divides this budget approximately equally
        # between both sides of the variant. Leave enough room for the full
        # longest prefix so this test exercises assembly rather than clipping.
        preferred_sequence_length=len(cdna_prefixes[-1]) * 2 + 1,
        variant_sequence_assembly=variant_sequence_assembly,
        min_assembly_overlap_size=1)
    variant_sequences = sequence_creator.reads_to_variant_sequences(
        variant=variant,
        reads=reads)
    protein_creator = ProteinSequenceCreator(
        protein_sequence_length=8,
        min_variant_sequence_coverage=2,
        min_transcript_prefix_length=3,
        max_transcript_mismatches=2)
    translations = protein_creator.all_pairs_translations(
        variant_sequences=variant_sequences,
        reference_contexts=[reference_context])

    if strand == "+":
        observed_candidates = {sequence.prefix for sequence in variant_sequences}
        compatible_sequence = "AAA"
    else:
        observed_candidates = {sequence.suffix for sequence in variant_sequences}
        compatible_sequence = reverse_complement_dna("AAA")
    eq_(observed_candidates,
        set(genomic_prefixes if strand == "+" else genomic_suffixes))
    eq_(len(translations), 1)
    translated_sequence = translations[0].untrimmed_variant_sequence
    eq_(translated_sequence.prefix if strand == "+" else translated_sequence.suffix,
        compatible_sequence)


def test_protein_sequence_creator_translates_coding_deletion():
    """
    Pin deletion handling through assembly, reference matching, and translation.
    """
    variant = Variant("17", 7676589, "CTC", "", grch38)
    transcript = grch38.transcripts_by_name("TP53-001")[0]
    protein_sequence_length = 10
    context_size = (protein_sequence_length + 1) * 3
    reference_key = ReferenceCodingSequenceKey.from_variant_and_transcript(
        variant=variant,
        transcript=transcript,
        context_size=context_size)

    # AlleleRead sequences use the positive genomic strand, whereas TP53 is
    # transcribed from the negative strand.
    prefix = reverse_complement_dna(
        reference_key.sequence_after_variant_locus)
    suffix = reverse_complement_dna(
        reference_key.sequence_before_variant_locus)
    reads = [
        AlleleRead(
            prefix=prefix,
            allele="",
            suffix=suffix,
            name="deletion_%d" % i)
        for i in range(4)
    ]

    creator = ProteinSequenceCreator(
        protein_sequence_length=protein_sequence_length,
        min_variant_sequence_coverage=2)
    translations = creator.translate_variant_reads(
        variant=variant,
        variant_reads=reads,
        transcript_id_whitelist={transcript.id})

    eq_(len(translations), 1)
    translation = translations[0]
    eq_(translation.amino_acids, "MEPQSD")
    eq_(translation.contains_mutation, True)
    eq_(translation.mutation_start_idx, 1)
    eq_(translation.mutation_end_idx, 1)
    eq_(translation.untrimmed_variant_sequence.alt, "")


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
