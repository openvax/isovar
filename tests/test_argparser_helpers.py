from .common import eq_
from isovar.cli.rna_args import (
    make_rna_reads_arg_parser,
    read_collector_from_args,
)
from isovar.cli.reference_context_args import add_reference_context_args
from isovar.cli.protein_sequence_args import add_protein_sequence_args
from isovar.cli.variant_sequences_args import add_variant_sequence_args


def test_extend_parser():
    parser = make_rna_reads_arg_parser()
    extra_arg_fns = [
        add_reference_context_args,
        add_protein_sequence_args,
        add_variant_sequence_args,
    ]
    for fn in extra_arg_fns:
        fn(parser)
    args = parser.parse_args([
        "--vcf", "ABC.vcf",
        "--maf", "ZZZ.maf",
        "--variant", "chr1", "39000", "C", "G",
        "--bam", "xyz.bam"])
    eq_(args.vcf, ["ABC.vcf"])
    eq_(args.maf, ["ZZZ.maf"])
    eq_(args.variant, [["chr1", "39000", "C", "G"]])
    eq_(args.bam, "xyz.bam")


def test_use_soft_clipped_bases_default_off():
    parser = make_rna_reads_arg_parser()
    args = parser.parse_args(["--vcf", "x.vcf", "--bam", "x.bam"])
    rc = read_collector_from_args(args)
    eq_(rc.use_soft_clipped_bases, False)


def test_use_soft_clipped_bases_flag_on():
    parser = make_rna_reads_arg_parser()
    args = parser.parse_args([
        "--vcf", "x.vcf",
        "--bam", "x.bam",
        "--use-soft-clipped-bases"])
    rc = read_collector_from_args(args)
    eq_(rc.use_soft_clipped_bases, True)
