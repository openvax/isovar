from isovar.args import (
    extend_parser_with_somatic_vcf_args,
    extend_parser_with_rna_args,
    extend_parser_with_reference_context_args,
)
from argparse import ArgumentParser
from nose.tools import eq_

def test_extend_parser():
    parser = ArgumentParser()
    fns = [
        extend_parser_with_somatic_vcf_args,
        extend_parser_with_rna_args,
        extend_parser_with_reference_context_args,
    ]
    for fn in fns:
        parser = fn(parser)
    args = parser.parse_args([
        "--vcf", "ABC.vcf",
        "--bam", "xyz.bam"])
    eq_(args.vcf, "ABC.vcf")
    eq_(args.bam, "xyz.bam")


if __name__ == "__main__":
    test_extend_parser()
