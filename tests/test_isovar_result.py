from isovar import run_isovar
from isovar import ProteinSequence
from varcode import Variant
from testing_helpers import data_path

from nose.tools import eq_


def test_isovar_result_property_types():
    for result in run_isovar(
            variants=data_path("data/b16.f10/b16.vcf"),
            alignment_file=data_path("data/b16.f10/b16.combined.sorted.bam")):
        # variant
        assert type(result.variant) is Variant

        # counts of genes and transcripts from variant
        assert type(result.num_overlapping_genes) is int
        assert type(result.num_overlapping_coding_genes) is int
        assert type(result.num_overlapping_transcripts) is int
        assert type(result.num_overlapping_coding_transcripts) is int

        # protein sequence
        assert type(result.top_protein_sequence) in (type(None), ProteinSequence)

        # counts of genes and transcripts from protein sequences
        assert type(result.num_genes_from_protein_sequences) is int
        assert type(result.num_genes_from_top_protein_sequence) is int
        assert type(result.num_transcripts_from_protein_sequences) is int
        assert type(result.num_transcripts_from_top_protein_sequence) is int

        # read and fragment counts
        assert type(result.num_ref_reads) is int
        assert type(result.num_alt_reads) is int
        assert type(result.num_other_reads) is int
        assert type(result.num_ref_fragments) is int
        assert type(result.num_alt_fragments) is int
        assert type(result.num_other_fragments) is int

        # read and fragment fractions
        assert type(result.fraction_ref_reads) is float
        assert type(result.fraction_alt_reads) is float
        assert type(result.fraction_other_reads) is float
        assert type(result.fraction_ref_fragments) is float
        assert type(result.fraction_alt_fragments) is float
        assert type(result.fraction_other_fragments) is float

        # read and fragment count ratios
        assert type(result.ratio_alt_to_other_reads) is float
        assert type(result.ratio_alt_to_other_fragments) is float
        assert type(result.ratio_other_to_alt_fragments) is float
        assert type(result.ratio_other_to_alt_reads) is float
        assert type(result.ratio_ref_to_other_fragments) is float
        assert type(result.ratio_other_to_ref_fragments) is float
        assert type(result.ratio_other_to_ref_reads) is float

        # this property aggregates all filters
        assert result.passes_all_filters in {True, False}

        assert type(result.protein_sequence_mutation_start_idx) in (int, type(None))
        assert type(result.protein_sequence_mutation_end_idx) in (int, type(None))


def test_isovar_result_nonsyn_variants():
    for result in run_isovar(
            variants=data_path("data/b16.f10/b16.vcf"),
            alignment_file=data_path("data/b16.f10/b16.combined.sorted.bam")):
        print(result.variant)
        print(result.predicted_effect)
        if result.has_mutant_protein_sequence_from_rna:
            assert result.num_amino_acid_mismatches_from_predicted_effect is not None
            assert result.num_amino_acid_mismatches_from_reference is not None
            assert result.num_amino_acid_mismatches_from_reference > 0
            eq_(result.protein_sequence_matches_reference, False)
            eq_(result.protein_sequence_contains_mutation, True)
        else:
            assert result.num_amino_acid_mismatches_from_predicted_effect is None
            assert result.num_amino_acid_mismatches_from_reference is None
            assert result.protein_sequence_matches_predicted_mutation_effect is None
            assert result.protein_sequence_matches_reference is None
            assert result.protein_sequence_contains_mutation is None


def test_isovar_result_clone():
    for result in run_isovar(
            variants=data_path("data/b16.f10/b16.vcf"),
            alignment_file=data_path("data/b16.f10/b16.combined.sorted.bam")):
        result2 = result.clone()
        eq_(result, result2)

def test_isovar_result_clone_with_updates():
    for result in run_isovar(
            variants=data_path("data/b16.f10/b16.vcf"),
            alignment_file=data_path("data/b16.f10/b16.combined.sorted.bam")):
        result2 = result.clone_with_updates(variant=None)
        assert result != result2


def test_isovar_result_str():
    for result in run_isovar(
            variants=data_path("data/b16.f10/b16.vcf"),
            alignment_file=data_path("data/b16.f10/b16.combined.sorted.bam")):
        s = str(result)
        assert len(s) > 0
        assert s.startswith("IsovarResult(")
        assert s.endswith(")")