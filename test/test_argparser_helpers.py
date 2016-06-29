from argparse import ArgumentParser
from isovar.args.variants import add_somatic_vcf_args
from isovar.args.rna_reads import add_rna_args
from isovar.args.reference_context import add_reference_context_args
from isovar.args.protein_sequences import add_protein_sequence_args
from isovar.args.variant_sequecnes import add_variant_sequence_args\

from isovar.protein_sequences import protein_sequences_dataframe_from_args

from testing_helpers import data_path

def test_variants_to_protein_sequences_dataframe_from_args():
    vcf_path = data_path("data/b16.f10/b16.vcf")
    bam_path = data_path("data/b16.f10/b16.combined.bam")

    parser = ArgumentParser()
    add_protein_sequence_args(parser)
    add_reference_context_args(parser)
    add_rna_args(parser)
    add_somatic_vcf_args(parser)
    add_variant_sequence_args(parser)
    args = parser.parse_args([
        "--vcf", vcf_path,
        "--bam", bam_path,
        "--genome", "mm10",
    ])
    df = protein_sequences_dataframe_from_args(args)
    assert len(df) > 0
