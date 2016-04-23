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

from nose.tools import eq_
from isovar.translation import Translation
from isovar.protein_sequence import (
    ProteinSequence,
    sort_protein_sequences,
    variants_to_protein_sequences_dataframe
)

from testing_helpers import load_bam, load_vcf


# fields of a ProteinSequence:
#   translations
#   supporting_variant_reads
#   total_variant_reads
#   supporting_transcripts
#   total_transcripts
#   gene

def make_dummy_translation(
        amino_acids="MKHW",  # ATG=M|AAA=K|CAC=H|TGG=W
        cdna_sequence="CCCATGAAACACTGGTAG",
        offset_to_first_complete_codon=3,
        variant_cdna_interval_start=8,  # assuming variant was AAC>AAA
        variant_cdna_interval_end=9,
        variant_aa_interval_start=1,
        variant_aa_interval_end=2,
        number_mismatches=1):
    return Translation(
        variant_sequence_in_reading_frame=None,
        cdna_sequence=cdna_sequence,
        offset_to_first_complete_codon=offset_to_first_complete_codon,
        variant_cdna_interval_start=variant_cdna_interval_start,
        variant_cdna_interval_end=variant_cdna_interval_end,
        reference_cdna_sequence_before_variant=cdna_sequence[:variant_cdna_interval_start],
        number_mismatches=number_mismatches,
        amino_acids=amino_acids,
        variant_aa_interval_start=variant_aa_interval_start,
        variant_aa_interval_end=variant_aa_interval_end,
        frameshift=False,
        ends_with_stop_codon=False,
        variant_sequence=None,
        reference_context=None)

def make_dummy_protein_sequence(
        n_supporting_variant_reads,
        n_supporting_variant_sequences,
        n_supporting_reference_transcripts,
        n_total_variant_sequences=None,
        n_total_variant_reads=None,
        n_total_reference_transcripts=None,
        gene=["TP53"],
        amino_acids="MKHW",  # ATG=M|AAA=K|CAC=H|TGG=W
        cdna_sequence="CCCATGAAACACTGGTAG",
        variant_cdna_interval_start=8,  # assuming variant was AAC>AAA
        variant_cdna_interval_end=9,
        variant_aa_interval_start=1,
        variant_aa_interval_end=2,
        number_mismatches=1):
    """
    Creates ProteinSequence object with None filled in for most fields
    """
    if n_total_variant_reads is None:
        n_total_variant_reads = n_supporting_variant_reads

    if n_total_variant_sequences is None:
        n_total_variant_sequences = n_supporting_variant_sequences

    if n_total_reference_transcripts is None:
        n_total_reference_transcripts = n_total_reference_transcripts

    assert n_supporting_variant_sequences <= n_supporting_variant_reads
    assert n_supporting_variant_sequences <= n_total_variant_sequences
    assert n_supporting_reference_transcripts <= n_total_reference_transcripts

    n_translations = n_total_reference_transcripts * n_total_variant_sequences

    translation = make_dummy_translation()

    return ProteinSequence(
        translations=[translation] * n_translations,
        supporting_variant_reads=[None] * n_supporting_variant_reads,
        total_variant_reads=n_total_variant_reads,
        supporting_transcripts=[None] * n_supporting_reference_transcripts,
        total_transcripts=n_total_reference_transcripts,
        gene=gene,
        amino_acids=amino_acids,
        variant_aa_interval_start=variant_aa_interval_start,
        variant_aa_interval_end=variant_aa_interval_end,
        ends_with_stop_codon=translation.ends_with_stop_codon,
        frameshift=translation.frameshift)

def test_sort_protein_sequences():
    protseq_most_reads = make_dummy_protein_sequence(
        n_supporting_variant_reads=50,
        n_supporting_variant_sequences=1,
        n_supporting_reference_transcripts=2,
        n_total_variant_sequences=3,
        n_total_variant_reads=100,
        n_total_reference_transcripts=5)

    protseq_most_reference_transcripts = make_dummy_protein_sequence(
        n_supporting_variant_reads=40,
        n_supporting_variant_sequences=1,
        n_supporting_reference_transcripts=3,
        n_total_variant_sequences=3,
        n_total_variant_reads=100,
        n_total_reference_transcripts=5)

    protseq_fewest_reads_or_transcripts = make_dummy_protein_sequence(
        n_supporting_variant_reads=10,
        n_supporting_variant_sequences=1,
        n_supporting_reference_transcripts=1,
        n_total_variant_sequences=3,
        n_total_variant_reads=100,
        n_total_reference_transcripts=5)
    unsorted_protein_sequences = [
        protseq_fewest_reads_or_transcripts,
        protseq_most_reads,
        protseq_most_reference_transcripts
    ]
    expected_order = [
        protseq_most_reads,
        protseq_most_reference_transcripts,
        protseq_fewest_reads_or_transcripts,
    ]
    eq_(sort_protein_sequences(unsorted_protein_sequences), expected_order)

"""
protein_sequence_length=PROTEIN_SEQUENCE_LEGNTH,
        min_reads_supporting_rna_sequence=MIN_READS_SUPPORTING_VARIANT_CDNA_SEQUENCE,
        min_transcript_prefix_length=MIN_TRANSCRIPT_PREFIX_LENGTH,
        max_transcript_mismatches=MAX_REFERENCE_TRANSCRIPT_MISMATCHES,
        max_protein_sequences_per_variant=MAX_PROTEIN_SEQUENCES_PER_VARIANT,
        min_mapping_quality=MIN_READ_MAPPING_QUALITY
"""

def test_variants_to_protein_sequences_dataframe_one_sequence_per_variant():
    variants = load_vcf("data/b16.f10/b16.vcf")
    samfile = load_bam("data/b16.f10/b16.combined.sorted.bam")
    df = variants_to_protein_sequences_dataframe(
        variants,
        samfile,
        min_mapping_quality=0,
        max_protein_sequences_per_variant=1)
    print(df)
    eq_(
        len(df),
        len(variants),
        "Expected %d entries, got %d" % (len(variants), len(df)))

def test_variants_to_protein_sequences_dataframe_filtered_all_reads_by_mapping_quality():
    # since the B16 BAM has all MAPQ=255 values then all the reads should get dropped
    # if we set the minimum quality to 256
    variants = load_vcf("data/b16.f10/b16.vcf")
    samfile = load_bam("data/b16.f10/b16.combined.sorted.bam")
    df = variants_to_protein_sequences_dataframe(
        variants,
        samfile,
        min_mapping_quality=256,
        max_protein_sequences_per_variant=1)
    print(df)
    eq_(
        len(df),
        0,
        "Expected 0 entries, got %d" % (len(df),))

def test_variants_to_protein_sequences_dataframe_protein_sequence_length():
    # make sure that the protein sequencces we get back are the length we asked for 
    variants = load_vcf("data/b16.f10/b16.vcf")
    samfile = load_bam("data/b16.f10/b16.combined.sorted.bam")

    for desired_length in range(9, 20, 3):
        df = variants_to_protein_sequences_dataframe(
            variants,
            samfile,
            max_protein_sequences_per_variant=1,
            protein_sequence_length=desired_length)
        eq_(
            len(df),
            len(variants),
            "Expected %d entries for protein_sequnece_length=%d, got %d results" % (
                len(variants),
                desired_length,
                len(df)))
        protein_sequences = df["amino_acids"]
        print(protein_sequences)
        protien_sequence_lengths = protein_sequences.str.len()
        assert (protien_sequence_lengths == desired_length).all(), protien_sequence_lengths

if __name__ == "__main__":
    test_sort_protein_sequences()
    test_variants_to_protein_sequences_dataframe_one_sequence_per_variant()
    test_variants_to_protein_sequences_dataframe_filtered_all_reads_by_mapping_quality()
    test_variants_to_protein_sequences_dataframe_protein_sequence_length()
