from __future__ import print_function, division, absolute_import

from nose.tools import eq_
from varcode import Variant

from isovar.allele_read import AlleleRead
from isovar.read_collector import ReadCollector

from mock_objects import MockAlignmentFile, make_pysam_read
from genomes_for_testing import grch38


def test_partitioned_read_sequences_snv():
    """
    test_partitioned_read_sequences_snv : Test that read gets correctly
    partitioned for chr1:4 T>G where the sequence for chr1 is assumed
    to be "ACCTTG"
    """
    # chr1_seq = "ACCTTG"
    chromosome = "1"
    location = 4
    ref = "T"
    alt = "G"

    variant = Variant(
        chromosome,
        location,
        ref,
        alt,
        grch38,
        normalize_contig_names=False)

    read = make_pysam_read(seq="ACCGTG", cigar="6M", mdtag="3G2")

    samfile = MockAlignmentFile(
        references=(chromosome,),
        reads=[read])
    read_creator = ReadCollector()
    variant_reads = read_creator.allele_reads_supporting_variant(
        variant=variant,
        alignment_file=samfile)
    print(variant_reads)
    assert len(variant_reads) == 1
    variant_read = variant_reads[0]
    expected = AlleleRead(
        name=read.qname,
        prefix="ACC",
        allele="G",
        suffix="TG")
    eq_(variant_read, expected)


def test_partitioned_read_sequences_insertion():
    """
    test_partitioned_read_sequences_insertion : Test that read gets correctly
    partitioned for chr1:4 T>TG
    where the sequence for chr1 is assumed to be "ACCTTG"
    and the variant sequence is "ACCTGTG"
    """
    # chr1_seq = "ACCTTG"
    chromosome = "1"
    location = 4
    ref = "T"
    alt = "TG"
    variant = Variant(
        chromosome, location, ref, alt, grch38, normalize_contig_names=False)

    read = make_pysam_read(
        seq=b"ACCTGTG",
        cigar="4M1I2M",
        mdtag="6")

    samfile = MockAlignmentFile(
        references=(chromosome,),
        reads=[read])
    read_creator = ReadCollector()

    variant_reads = read_creator.allele_reads_supporting_variant(
        alignment_file=samfile,
        variant=variant)
    print(variant_reads)
    assert len(variant_reads) == 1
    variant_read = variant_reads[0]
    expected = AlleleRead(
        name=read.qname,
        prefix="ACCT",
        allele="G",
        suffix="TG")
    eq_(variant_read, expected)


def test_partitioned_read_sequences_deletion():
    """
    test_partitioned_read_sequences_deletion : Test that read gets correctly
    partitioned for chr1:4 TT>T where the sequence for chr1 is assumed to
    be "ACCTTG"
    """
    # chr1_seq = "ACCTTG"
    chromosome = "1"
    location = 4
    ref = "TT"
    alt = "T"
    variant = Variant(
        chromosome, location, ref, alt, grch38, normalize_contig_names=False)

    read = make_pysam_read(
        seq="ACCTG",
        cigar="4M1D1M",
        mdtag="4^T1")
    samfile = MockAlignmentFile(
        references=(chromosome,),
        reads=[read])
    read_creator = ReadCollector()
    variant_reads = read_creator.allele_reads_supporting_variant(
        alignment_file=samfile,
        variant=variant)
    print(variant_reads)
    assert len(variant_reads) == 1
    variant_read = variant_reads[0]
    expected = AlleleRead(
        name=read.qname,
        prefix="ACCT",
        allele="",
        suffix="G")
    eq_(variant_read, expected)
