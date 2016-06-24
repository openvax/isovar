# Copyright (c) 2016. Mount Sinai School of Medicine
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

from __future__ import print_function, division, absolute_import

from nose.tools import eq_
from varcode import Variant
from isovar.locus_reads import (
    LocusRead,
    locus_read_generator,
    locus_reads_dataframe,
)

from mock_read_data import DummySamFile, make_read
from testing_helpers import assert_equal_fields, load_bam, data_path

def test_locus_reads_snv():
    """
    test_partitioned_read_sequences_snv : Test that read gets correctly
    partitioned for chr1:4 T>G where the sequence for chr1 is assumed
    to be "ACCTTG"
    """
    # chr1_seq = "ACCTTG"
    variant = Variant(
        "chromosome",
        4,
        ref="T",
        alt="G",
        normalize_contig_name=False)

    pysam_read = make_read(seq="ACCGTG", cigar="6M", mdtag="3G2")

    samfile = DummySamFile(reads=[pysam_read])
    reads = list(locus_read_generator(
        samfile=samfile,
        chromosome="chromosome",
        base1_position_before_variant=variant.start - 1,
        base1_position_after_variant=variant.start + 1))
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
        base0_read_position_before_variant=2,
        base0_read_position_after_variant=4)
    assert_equal_fields(read, expected)

def test_locus_reads_insertion():
    """
    test_partitioned_read_sequences_insertion : Test that read gets correctly
    partitioned for chr1:4 T>TG
    where the sequence for chr1 is assumed to be "ACCTTG"
    and the variant sequence is "ACCTGTG"
    """
    variant = Variant(
        "chromosome", 4, ref="T", alt="TG", normalize_contig_name=False)

    pysam_read = make_read(seq="ACCTGTG", cigar="4M1I2M", mdtag="6")

    samfile = DummySamFile(reads=[pysam_read])
    reads = list(locus_read_generator(
        samfile=samfile,
        chromosome="chromosome",
        base1_position_before_variant=variant.start,
        base1_position_after_variant=variant.start + 1))
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
        base0_read_position_before_variant=3,
        base0_read_position_after_variant=5)
    assert_equal_fields(read, expected)

def test_locus_reads_deletion():
    """
    test_partitioned_read_sequences_deletion : Test that read gets correctly
    partitioned for chr1:4 TT>T where the sequence for chr1 is assumed to
    be "ACCTTG"
    """
    # normalization of this variant will turn it into the deletion of
    # "T" at base-1 position 5
    variant = Variant(
        "chromosome", 4, ref="TT", alt="T", normalize_contig_name=False)
    print(variant)
    pysam_read = make_read(seq="ACCTG", cigar="4M1D1M", mdtag="4^T1")

    samfile = DummySamFile(reads=[pysam_read])
    reads = list(locus_read_generator(
        samfile=samfile,
        chromosome="chromosome",
        base1_position_before_variant=variant.start - 1,
        base1_position_after_variant=variant.start + 1))
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
        base0_read_position_before_variant=3,
        base0_read_position_after_variant=4)
    assert_equal_fields(read, expected)

def test_locus_reads_substitution_longer():
    # test C>GG subsitution at second nucleotide of reference sequence "ACCTTG",
    # the alignment is interpreted as a C>G variant followed by an insertion of
    # another G
    variant = Variant(
        "chromosome", 2, ref="C", alt="GG", normalize_contig_name=False)
    print(variant)
    pysam_read = make_read(seq="AGGCTTG", cigar="2M1I4M", mdtag="1C4")

    samfile = DummySamFile(reads=[pysam_read])
    reads = list(locus_read_generator(
        samfile=samfile,
        chromosome="chromosome",
        base1_position_before_variant=1,
        base1_position_after_variant=3))
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
        base0_read_position_before_variant=0,
        base0_read_position_after_variant=3)
    assert_equal_fields(read, expected)

def test_locus_reads_substitution_shorter():
    # test CC>G subsitution at 2nd and 3rd nucleotides of reference sequence
    # "ACCTTG", for which the alignment is interpreted as a C>G variant
    # followed by the deletion of a C
    variant = Variant(
        "chromosome", 2, ref="CC", alt="G", normalize_contig_name=False)
    print(variant)
    pysam_read = make_read(seq="AGTTG", cigar="2M1D3M", mdtag="1C^C4")

    samfile = DummySamFile(reads=[pysam_read])
    reads = list(locus_read_generator(
        samfile=samfile,
        chromosome="chromosome",
        base1_position_before_variant=1,
        base1_position_after_variant=4))
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
        base0_read_position_before_variant=0,
        base0_read_position_after_variant=2)
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
        samfile=sam_all_variants,
        chromosome="chr4",
        base1_position_before_variant=45802538,
        base1_position_after_variant=45802540)
    print(df)
    eq_(len(df), n_reads_expected)
