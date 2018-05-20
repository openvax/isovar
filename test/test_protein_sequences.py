# Copyright (c) 2016-2018. Mount Sinai School of Medicine
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

from isovar.cli.protein_sequence_args import (
    protein_sequences_dataframe_from_args,
    make_protein_sequences_arg_parser,
)
from isovar.protein_sequences import (
    ProteinSequence,
    sort_protein_sequences,
    reads_generator_to_protein_sequences_generator,
    protein_sequences_generator_to_dataframe,
)
from isovar.allele_reads import reads_overlapping_variants
from isovar.variant_sequence_in_reading_frame import VariantSequenceInReadingFrame
from varcode import VariantCollection

from testing_helpers import load_bam, load_vcf, data_path


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
    varseq_in_orf = VariantSequenceInReadingFrame(
        cdna_sequence=cdna_sequence,
        offset_to_first_complete_codon=offset_to_first_complete_codon,
        variant_cdna_interval_start=variant_cdna_interval_start,
        variant_cdna_interval_end=variant_cdna_interval_end,
        reference_cdna_sequence_before_variant=cdna_sequence[:variant_cdna_interval_start],
        reference_cdna_sequence_after_variant=cdna_sequence[variant_cdna_interval_end:],
        number_mismatches_before_variant=number_mismatches,
        number_mismatches_after_variant=0)
    return Translation(
        variant_sequence_in_reading_frame=varseq_in_orf,
        amino_acids=amino_acids,
        variant_aa_interval_start=variant_aa_interval_start,
        variant_aa_interval_end=variant_aa_interval_end,
        frameshift=False,
        ends_with_stop_codon=False,
        untrimmed_variant_sequence=None,
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
        overlapping_reads=[None] * n_total_variant_reads,
        ref_reads=[],
        alt_reads=[None] * n_total_variant_reads,
        alt_reads_supporting_protein_sequence=[None] * n_supporting_variant_reads,
        transcripts_supporting_protein_sequence=[None] * n_supporting_reference_transcripts,
        transcripts_overlapping_variant=[None] * n_supporting_reference_transcripts,
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


def variants_to_protein_sequences_dataframe(
        expressed_vcf="data/b16.f10/b16.expressed.vcf",
        not_expressed_vcf="data/b16.f10/b16.not-expressed.vcf",
        tumor_rna_bam="data/b16.f10/b16.combined.sorted.bam",
        min_mapping_quality=0,
        max_protein_sequences_per_variant=1,
        variant_sequence_assembly=False):
    """
    Helper function to load pair of VCFs and tumor RNA BAM
    and use them to generate a DataFrame of expressed variant protein
    sequences.
    """
    expressed_variants = load_vcf(expressed_vcf)
    not_expressed_variants = load_vcf(not_expressed_vcf)

    combined_variants = VariantCollection(
        list(expressed_variants) + list(not_expressed_variants))
    samfile = load_bam(tumor_rna_bam)

    allele_reads_generator = reads_overlapping_variants(
        variants=combined_variants,
        samfile=samfile,
        min_mapping_quality=min_mapping_quality)

    protein_sequences_generator = reads_generator_to_protein_sequences_generator(
        allele_reads_generator,
        max_protein_sequences_per_variant=max_protein_sequences_per_variant,
        variant_sequence_assembly=variant_sequence_assembly)
    df = protein_sequences_generator_to_dataframe(protein_sequences_generator)
    return df, expressed_variants, combined_variants


def test_variants_to_protein_sequences_dataframe_one_sequence_per_variant_with_assembly():
    df, expressed_variants, combined_variants = \
        variants_to_protein_sequences_dataframe(variant_sequence_assembly=True)
    print(df)
    eq_(len(df),
        len(expressed_variants),
        "Expected %d/%d entries to have RNA support, got %d" % (
            len(expressed_variants),
            len(combined_variants),
            len(df)))


def test_variants_to_protein_sequences_dataframe_one_sequence_per_variant_without_assembly():
    df, expressed_variants, combined_variants = \
        variants_to_protein_sequences_dataframe(variant_sequence_assembly=False)
    print(df)
    eq_(len(df),
        len(expressed_variants),
        "Expected %d/%d entries to have RNA support, got %d" % (
            len(expressed_variants),
            len(combined_variants),
            len(df)))


def test_variants_to_protein_sequences_dataframe_filtered_all_reads_by_mapping_quality():
    # since the B16 BAM has all MAPQ=255 values then all the reads should get dropped
    # if we set the minimum quality to 256
    variants = load_vcf("data/b16.f10/b16.vcf")
    samfile = load_bam("data/b16.f10/b16.combined.sorted.bam")
    allele_reads_generator = reads_overlapping_variants(
        variants=variants,
        samfile=samfile,
        min_mapping_quality=256)
    protein_sequences_generator = reads_generator_to_protein_sequences_generator(
        allele_reads_generator,
        max_protein_sequences_per_variant=1)
    df = protein_sequences_generator_to_dataframe(protein_sequences_generator)
    print(df)
    eq_(
        len(df),
        0,
        "Expected 0 entries, got %d: %s" % (len(df), df))


def test_variants_to_protein_sequences_dataframe_protein_sequence_length():
    expressed_variants = load_vcf("data/b16.f10/b16.expressed.vcf")
    parser = make_protein_sequences_arg_parser()
    parser.print_help()
    for desired_length in range(9, 20, 3):
        args = parser.parse_args([
            "--vcf", data_path("data/b16.f10/b16.vcf"),
            "--bam", data_path("data/b16.f10/b16.combined.sorted.bam"),
            "--max-protein-sequences-per-variant", "1",
            "--protein-sequence-length", str(desired_length),
        ])
        df = protein_sequences_dataframe_from_args(args)
        eq_(
            len(df),
            len(expressed_variants),
            "Expected %d entries for protein_sequence_length=%d, got %d results: %s" % (
                len(expressed_variants),
                desired_length,
                len(df),
                df))
        protein_sequences = df["amino_acids"]
        print(protein_sequences)
        protein_sequence_lengths = protein_sequences.str.len()
        assert (protein_sequence_lengths == desired_length).all(), (
            protein_sequence_lengths,)
