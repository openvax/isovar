from __future__ import print_function, division, absolute_import

from pyensembl import ensembl_grch38

from varcode import Variant
from nose.tools import eq_

from testing_helpers import load_bam

from isovar.variant_read import reads_supporting_variant


def test_partition_variant_reads_snv():
    samfile = load_bam("data/cancer-wgs-primary.chr12.bam")
    chromosome = "chr12"
    base1_location = 65857041
    ref = "G"
    alt = "C"
    variant = Variant(
        contig=chromosome,
        start=base1_location,
        ref=ref,
        alt=alt,
        ensembl=ensembl_grch38)
    variant_reads = reads_supporting_variant(
        samfile=samfile,
        chromosome=chromosome,
        variant=variant)
    assert len(variant_reads) > 1
    for variant_read in variant_reads:
        eq_(variant_read.allele, alt)


def test_partition_variant_reads_deletion():
    samfile = load_bam("data/cancer-wgs-primary.chr12.bam")
    chromosome = "chr12"
    base1_location = 70091490
    ref = "TTGTAGATGCTGCCTCTCC"
    alt = ""
    variant = Variant(
        contig=chromosome,
        start=base1_location,
        ref=ref,
        alt=alt,
        ensembl=ensembl_grch38)
    variant_reads = reads_supporting_variant(
        samfile=samfile,
        chromosome=chromosome,
        variant=variant)
    assert len(variant_reads) > 1
    for variant_read in variant_reads:
        eq_(variant_read.allele, alt)

if __name__ == "__main__":
    test_partition_variant_reads_snv()
    test_partition_variant_reads_deletion()
