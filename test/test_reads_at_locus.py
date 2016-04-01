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
from isovar.reads_at_locus import ReadAtLocus, gather_reads_at_locus

from mock_read_data import DummySamFile, make_read


def test_reads_at_locus_snv():
    """
    test_partitioned_read_sequences_snv : Test that read gets correctly
    partitioned for chr1:4 T>G where the sequence for chr1 is assumed
    to be "ACCTTG"
    """
    # chr1_seq = "ACCTTG"
    variant = Variant(
        "chromosome", 4, ref="T", alt="G", normalize_contig_name=False)

    pysam_read = make_read(seq="ACCGTG", cigar="6M", mdtag="3G2")

    samfile = DummySamFile(reads=[pysam_read])
    reads = list(gather_reads_at_locus(
        samfile=samfile,
        chromosome="chromosome",
        base1_position_before_variant=variant.start - 1,
        base1_position_after_variant=variant.start + 1))
    eq_(len(reads), 1)
    print(reads)
    read = reads[0]
    expected = ReadAtLocus(
        name=pysam_read.qname,
        sequence="ACCGTG",
        reference_positions=[0, 1, 2, 3, 4, 5],
        base_qualities=pysam_read.query_qualities,
        offset_before_variant=2,
        offset_after_variant=4)
    eq_(read, expected)

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
    eq_(len(reads), 1)
    print(reads)
    read = reads[0]
    expected = ReadAtLocus(
        name=pysam_read.qname,
        sequence="ACCTGTG",
        reference_positions=[0, 1, 2, 3, None, 4, 5],
        base_qualities=pysam_read.query_qualities,
        offset_before_variant=3,
        offset_after_variant=4)
    eq_(read, expected)

def test_reads_at_locus_deletion():
    """
    test_partitioned_read_sequences_deletion : Test that read gets correctly
    partitioned for chr1:4 TT>T where the sequence for chr1 is assumed to
    be "ACCTTG"
    """
    variant = Variant(
        "chromosome", 4, ref="TT", alt="T", normalize_contig_name=False)

    pysam_read = make_read(seq="ACCTG", cigar="4M1D1M", mdtag="4^T1")

    samfile = DummySamFile(reads=[pysam_read])
    reads = list(gather_reads_at_locus(
        samfile=samfile,
        chromosome="chromosome",
        base1_position_before_variant=variant.start - 1,
        base1_position_after_variant=variant.start + 1))
    eq_(len(reads), 1)
    print(reads)
    read = reads[0]
    expected = ReadAtLocus(
        name=pysam_read.qname,
        sequence="ACCTG",
        reference_positions=[0, 1, 2, 3, 5],
        base_qualities=pysam_read.query_qualities,
        offset_before_variant=2,
        offset_after_variant=4)
    eq_(read, expected)

def test_reads_at_locus_substitution_longer():
    # test C>TT subsitution
    pass

def test_reads_at_locus_substitution_shorter():
    # test TT>C substitution
    pass
