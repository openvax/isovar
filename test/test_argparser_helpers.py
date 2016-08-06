from nose.tools import eq_
from isovar.cli.rna_reads import make_rna_reads_arg_parser
from isovar.cli.reference_context import add_reference_context_args
from isovar.cli.protein_sequences import add_protein_sequence_args
from isovar.cli.variant_sequences import add_variant_sequence_args


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
