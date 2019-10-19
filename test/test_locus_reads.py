from __future__ import print_function, division, absolute_import

from nose.tools import eq_
from varcode import Variant
from isovar.locus_read import LocusRead
from isovar.dataframe_helpers import locus_reads_dataframe
from isovar.read_collector import ReadCollector

from mock_objects import MockAlignmentFile, make_pysam_read
from testing_helpers import assert_equal_fields, load_bam, data_path


def test_locus_reads_snv():
    """
    test_partitioned_read_sequences_snv : Test that read gets correctly
    partitioned for chr1:4 T>G where the sequence for chr1 is assumed
    to be "ACCTTG"
    """
    # chr1_seq = "ACCTTG"
    variant = Variant(
        "1",
        4,
        ref="T",
        alt="G")

    pysam_read = make_pysam_read(
        seq="ACCGTG",
        cigar="6M",
        mdtag="3G2")

    samfile = MockAlignmentFile(
        references=("chromosome",),
        reads=[pysam_read])
    read_creator = ReadCollector()
    reads = read_creator.get_locus_reads(samfile, "chromosome", variant.start - 1, variant.start)
    print(reads)
    assert len(reads) == 1, \
        "Expected to get back one read but instead got %d" % (
            len(reads),)
    read = reads[0]
    expected = LocusRead(
        name=pysam_read.qname,
        sequence=pysam_read.query_sequence,
        reference_positions=[0, 1, 2, 3, 4, 5],
        quality_scores=pysam_read.query_qualities,
        reference_base0_start_inclusive=3,
        reference_base0_end_exclusive=4,
        read_base0_start_inclusive=3,
        read_base0_end_exclusive=4)
    assert_equal_fields(read, expected)


def test_locus_reads_insertion():
    """
    test_partitioned_read_sequences_insertion : Test that read gets correctly
    partitioned for chr1:4 T>TG
    where the sequence for chr1 is assumed to be "ACCTTG"
    and the variant sequence is "ACCTGTG"
    """
    variant = Variant("1", 4, ref="T", alt="TG")

    pysam_read = make_pysam_read(seq="ACCTGTG", cigar="4M1I2M", mdtag="6")

    samfile = MockAlignmentFile(
        references={"chromosome"},
        reads=[pysam_read])
    read_creator = ReadCollector()
    reads = read_creator.get_locus_reads(
        samfile,
        "chromosome",
        variant.start,
        variant.start)
    print(reads)
    assert len(reads) == 1, \
        "Expected to get back one read but instead got %d" % (
            len(reads),)
    read = reads[0]
    expected = LocusRead(
        name=pysam_read.qname,
        sequence=pysam_read.query_sequence,
        # expect the inserted nucleotide to be missing a corresponding
        # ref position
        reference_positions=[0, 1, 2, 3, None, 4, 5],
        quality_scores=pysam_read.query_qualities,
        read_base0_start_inclusive=4,
        read_base0_end_exclusive=5,
        reference_base0_start_inclusive=4,
        reference_base0_end_exclusive=4)
    print("Actual: %s" % (read,))
    print("Expected: %s" % (expected,))
    assert_equal_fields(read, expected)


def test_locus_reads_deletion():
    """
    test_partitioned_read_sequences_deletion : Test that read gets correctly
    partitioned for chr1:4 TT>T where the sequence for chr1 is assumed to
    be "ACCTTG"
    """
    # normalization of this variant will turn it into the deletion of
    # "T" at base-1 position 5
    variant = Variant("1", 4, ref="TT", alt="T")
    pysam_read = make_pysam_read(seq="ACCTG", cigar="4M1D1M", mdtag="4^T1")

    samfile = MockAlignmentFile(
        references={"chromosome"},
        reads=[pysam_read])
    read_creator = ReadCollector()
    reads = read_creator.get_locus_reads(
        samfile,
        "chromosome",
        variant.start - 1,
        variant.start)
    print(reads)
    assert len(reads) == 1, \
        "Expected to get back one read but instead got %d" % (
            len(reads),)
    read = reads[0]
    expected = LocusRead(
        name=pysam_read.qname,
        sequence=pysam_read.query_sequence,
        reference_positions=[0, 1, 2, 3, 5],
        quality_scores=pysam_read.query_qualities,
        # missing would have gone after 4th nucleotide in the read
        read_base0_start_inclusive=4,
        read_base0_end_exclusive=4,
        reference_base0_start_inclusive=4,
        reference_base0_end_exclusive=5)
    assert_equal_fields(read, expected)


def test_locus_reads_substitution_longer():
    # test C>GG subsitution at second nucleotide of reference sequence "ACCTTG",
    # the alignment is interpreted as a C>G variant followed by an insertion of
    # another G
    variant = Variant("1", 2, ref="C", alt="GG")
    print(variant)
    pysam_read = make_pysam_read(seq="AGGCTTG", cigar="2M1I4M", mdtag="1C4")

    samfile = MockAlignmentFile(
        references={"chromosome"},
        reads=[pysam_read])
    read_creator = ReadCollector()
    reads = read_creator.get_locus_reads(
        samfile,
        "chromosome",
        1,
        2)
    print(reads)
    assert len(reads) == 1, \
        "Expected to get back one read but instead got %d" % (
            len(reads),)
    read = reads[0]
    expected = LocusRead(
        name=pysam_read.qname,
        sequence=pysam_read.query_sequence,
        reference_positions=[0, 1, None, 2, 3, 4, 5],
        quality_scores=pysam_read.query_qualities,
        read_base0_start_inclusive=1,
        read_base0_end_exclusive=3,
        reference_base0_start_inclusive=1,
        reference_base0_end_exclusive=2)
    assert_equal_fields(read, expected)


def test_locus_reads_substitution_shorter():
    # test CC>G subsitution at 2nd and 3rd nucleotides of reference sequence
    # "ACCTTG", for which the alignment is interpreted as a C>G variant
    # followed by the deletion of a C
    variant = Variant("1", 2, ref="CC", alt="G")
    print(variant)
    pysam_read = make_pysam_read(seq="AGTTG", cigar="2M1D3M", mdtag="1C^C4")

    samfile = MockAlignmentFile(
        references={"chromosome"},
        reads=[pysam_read])
    read_creator = ReadCollector()
    reads = read_creator.get_locus_reads(
        samfile,
        "chromosome",
        1,
        3)
    assert len(reads) == 1, \
        "Expected to get back one read but instead got %d" % (
            len(reads),)
    print(reads)
    read = reads[0]
    expected = LocusRead(
        name=pysam_read.qname,
        sequence=pysam_read.query_sequence,
        reference_positions=[0, 1, 3, 4, 5],
        quality_scores=pysam_read.query_qualities,
        read_base0_start_inclusive=1,
        read_base0_end_exclusive=2,
        reference_base0_start_inclusive=1,
        reference_base0_end_exclusive=3)
    assert_equal_fields(read, expected)


def test_locus_reads_dataframe():
    sam_all_variants = load_bam("data/b16.f10/b16.combined.bam")

    n_reads_expected = 0

    sam_path_single_variant = data_path(
        "data/b16.f10/b16.f10.127a.aldh1b1.chr4.45802539.refG.altC.sam")
    with open(sam_path_single_variant) as f:
        for line in f:
            if line.startswith("HWI"):
                n_reads_expected += 1
    # we know from inspecting the file that *one* of the reads overlapping this
    # variant has a CIGAR string of N at the location before and thus we'll
    # be missing that read.
    #
    # TODO: figure out what to do when the variant nucleotide is at the start or
    # end of an exon, since that won't have mapping positions on both its left
    # and right
    n_reads_expected -= 1

    print("Found %d sequences in %s" % (n_reads_expected, sam_path_single_variant))
    df = locus_reads_dataframe(
        alignments=sam_all_variants,
        chromosome="chr4",
        base0_start=45802538,
        base0_end=45802539)
    print(df)
    eq_(len(df), n_reads_expected)
