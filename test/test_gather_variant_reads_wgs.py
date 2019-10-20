from __future__ import print_function, division, absolute_import

from pyensembl import ensembl_grch38

from varcode import Variant
from nose.tools import eq_

from testing_helpers import load_bam

from isovar import ReadCollector


def test_partition_variant_reads_snv():
    alignment_file = load_bam("data/cancer-wgs-primary.chr12.bam")
    chromosome = "12"
    base1_location = 65857041
    ref = "G"
    alt = "C"
    variant = Variant(
        contig=chromosome,
        start=base1_location,
        ref=ref,
        alt=alt,
        ensembl=ensembl_grch38)
    read_collector = ReadCollector()
    read_evidence = read_collector.read_evidence_for_variant(
        alignment_file=alignment_file,
        variant=variant)
    alt_reads = read_evidence.alt_reads
    assert len(alt_reads) > 1
    for variant_read in alt_reads:
        eq_(variant_read.allele, alt)


def test_partition_variant_reads_deletion():
    alignment_file = load_bam("data/cancer-wgs-primary.chr12.bam")
    chromosome = "12"
    base1_location = 70091490
    ref = "TTGTAGATGCTGCCTCTCC"
    alt = ""
    variant = Variant(
        contig=chromosome,
        start=base1_location,
        ref=ref,
        alt=alt,
        ensembl=ensembl_grch38)
    read_collector = ReadCollector()
    read_evidence = read_collector.read_evidence_for_variant(
        alignment_file=alignment_file,
        variant=variant)
    assert len(read_evidence.alt_reads) > 1
    for variant_read in read_evidence.alt_reads:
        eq_(variant_read.allele, alt)

if __name__ == "__main__":
    test_partition_variant_reads_snv()
    test_partition_variant_reads_deletion()
