from nose.tools import eq_
from argparse import ArgumentParser
from isovar.args.variants import add_somatic_vcf_args
from isovar.args.rna_reads import add_rna_args
from isovar.args.reference_context import add_reference_context_args
from isovar.args.protein_sequences import add_protein_sequence_args
from isovar.args.variant_sequences import add_variant_sequence_args


def test_extend_parser():
    parser = ArgumentParser()
    fns = [
        add_somatic_vcf_args,
        add_rna_args,
        add_reference_context_args,
        add_protein_sequence_args,
        add_variant_sequence_args,
    ]
    for fn in fns:
        fn(parser)
    args = parser.parse_args([
        "--vcf", "ABC.vcf",
        "--bam", "xyz.bam"])
    eq_(args.vcf, "ABC.vcf")
    eq_(args.bam, "xyz.bam")
