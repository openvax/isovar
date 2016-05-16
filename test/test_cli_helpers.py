from argparse import ArgumentParser
from isovar.args import (
    add_somatic_vcf_args,
    add_rna_args,
    add_reference_context_args,
    add_protein_sequence_args,
)
from isovar.cli_helpers import variants_to_protein_sequences_dataframe_from_args

from testing_helpers import data_path

def test_variants_to_protein_sequences_dataframe_from_args():
    vcf_path = data_path("data/b16.f10/b16.vcf")
    bam_path = data_path("data/b16.f10/b16.combined.bam")

    parser = \
        add_protein_sequence_args(
            add_reference_context_args(
                add_rna_args(
                    add_somatic_vcf_args(ArgumentParser()))))
    args = parser.parse_args([
        "--vcf", vcf_path,
        "--bam", bam_path,
        "--genome", "mm10",
    ])
    df = variants_to_protein_sequences_dataframe_from_args(args)
    assert len(df) > 0
