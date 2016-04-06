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

from varcode import Variant
from isovar.reads_at_locus import ReadAtLocus, gather_reads_at_locus

from mock_read_data import DummySamFile, make_read
from testing_helpers import assert_equal_fields

def test_reads_at_locus_snv():
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
    reads = list(gather_reads_at_locus(
        samfile=samfile,
        chromosome="chromosome",
        base1_position_before_variant=variant.start - 1,
        base1_position_after_variant=variant.start + 1))
    print(reads)
    assert len(reads) == 1, \
        "Expected to get back one read but instead got %d" % (
            len(reads),)
    read = reads[0]
    expected = ReadAtLocus(
        name=pysam_read.qname,
        sequence=pysam_read.query_sequence.decode('ascii'),
        reference_positions=[0, 1, 2, 3, 4, 5],
        quality_scores=pysam_read.query_qualities,
        base0_read_position_before_variant=2,
        base0_read_position_after_variant=4)
    assert_equal_fields(read, expected)

def test_reads_at_locus_insertion():
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
    reads = list(gather_reads_at_locus(
        samfile=samfile,
        chromosome="chromosome",
        base1_position_before_variant=variant.start,
        base1_position_after_variant=variant.start + 1))
    print(reads)
    assert len(reads) == 1, \
        "Expected to get back one read but instead got %d" % (
            len(reads),)
    read = reads[0]
    expected = ReadAtLocus(
        name=pysam_read.qname,
        sequence=pysam_read.query_sequence.decode('ascii'),
        # expect the inserted nucleotide to be missing a corresponding
        # ref position
        reference_positions=[0, 1, 2, 3, None, 4, 5],
        quality_scores=pysam_read.query_qualities,
        base0_read_position_before_variant=3,
        base0_read_position_after_variant=5)
    assert_equal_fields(read, expected)

def test_reads_at_locus_deletion():
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
    reads = list(gather_reads_at_locus(
        samfile=samfile,
        chromosome="chromosome",
        base1_position_before_variant=variant.start - 1,
        base1_position_after_variant=variant.start + 1))
    print(reads)
    assert len(reads) == 1, \
        "Expected to get back one read but instead got %d" % (
            len(reads),)
    read = reads[0]
    expected = ReadAtLocus(
        name=pysam_read.qname,
        sequence=pysam_read.query_sequence.decode('ascii'),
        reference_positions=[0, 1, 2, 3, 5],
        quality_scores=pysam_read.query_qualities,
        base0_read_position_before_variant=3,
        base0_read_position_after_variant=4)
    assert_equal_fields(read, expected)

def test_reads_at_locus_substitution_longer():
    # test C>GG subsitution at second nucleotide of reference sequence "ACCTTG",
    # the alignment is interpreted as a C>G variant followed by an insertion of
    # another G
    variant = Variant(
        "chromosome", 2, ref="C", alt="GG", normalize_contig_name=False)
    print(variant)
    pysam_read = make_read(seq="AGGCTTG", cigar="2M1I4M", mdtag="1C4")

    samfile = DummySamFile(reads=[pysam_read])
    reads = list(gather_reads_at_locus(
        samfile=samfile,
        chromosome="chromosome",
        base1_position_before_variant=1,
        base1_position_after_variant=3))
    print(reads)
    assert len(reads) == 1, \
        "Expected to get back one read but instead got %d" % (
            len(reads),)
    read = reads[0]
    expected = ReadAtLocus(
        name=pysam_read.qname,
        sequence=pysam_read.query_sequence.decode('ascii'),
        reference_positions=[0, 1, None, 2, 3, 4, 5],
        quality_scores=pysam_read.query_qualities,
        base0_read_position_before_variant=0,
        base0_read_position_after_variant=3)
    assert_equal_fields(read, expected)

def test_reads_at_locus_substitution_shorter():
    # test CC>G subsitution at 2nd and 3rd nucleotides of reference sequence
    # "ACCTTG", for which the alignment is interpreted as a C>G variant
    # followed by the deletion of a C
    variant = Variant(
        "chromosome", 2, ref="CC", alt="G", normalize_contig_name=False)
    print(variant)
    pysam_read = make_read(seq="AGTTG", cigar="2M1D3M", mdtag="1C^C4")

    samfile = DummySamFile(reads=[pysam_read])
    reads = list(gather_reads_at_locus(
        samfile=samfile,
        chromosome="chromosome",
        base1_position_before_variant=1,
        base1_position_after_variant=4))
    assert len(reads) == 1, \
        "Expected to get back one read but instead got %d" % (
            len(reads),)
    print(reads)
    read = reads[0]
    expected = ReadAtLocus(
        name=pysam_read.qname,
        sequence=pysam_read.query_sequence.decode('ascii'),
        reference_positions=[0, 1, 3, 4, 5],
        quality_scores=pysam_read.query_qualities,
        base0_read_position_before_variant=0,
        base0_read_position_after_variant=2)
    assert_equal_fields(read, expected)
