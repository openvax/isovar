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
import pysam
from varcode import Variant
from isovar.variant_reads import gather_reads_for_single_variant

class DummyPileupElement(object):
    def __init__(self, alignment, is_refskip, is_del):
        self.alignment = alignment
        self.is_del = is_del
        self.is_refskip = is_refskip

class DummyPileupColumn(object):
    def __init__(self, pos, reads, is_del=False, is_refskip=False):
        self.pos = pos
        self.pileups = [
            DummyPileupElement(read, is_del=is_del, is_refskip=is_refskip)
            for read in reads]

class DummySamFile(object):
    """
    Used instead of real AlignmentFile objects for test.
    """
    def __init__(self, reads):
        self.reads = reads

    def fetch(self, *args, **kwargs):
        return self.reads

    def pileup(self, chromosome, start, end):
        for i in range(start, end + 1):
            yield DummyPileupColumn(pos=i + 1, reads=self.reads)

def make_read(seq, cigar, mdtag=None, name="dummy", mapq=10, baseq=30):
    read = pysam.AlignedSegment()
    read.seq = seq
    read.cigarstring = cigar
    if mdtag:
        read.set_tag("MD", mdtag)
    read.qname = name
    read.mapq = mapq
    read.qual = bytes(
        pysam.qualities_to_qualitystring([baseq] * len(seq)), 'ascii')
    return read

def test_partitioned_read_sequences_snv():
    """
    test_partitioned_read_sequences_snv : Test that read gets correctly
    partitioned for chr1:4 T>G where the sequence for chr1 is assumed
    to be "ACCTTG"
    """
    # chr1_seq = "ACCTTG"
    chromosome = "chromosome"
    location = 4
    ref = "T"
    alt = "G"

    variant = Variant(
        chromosome, location, ref, alt, normalize_contig_name=False)

    read = make_read(seq="ACCGTG", cigar="6M", mdtag="3G2")

    samfile = DummySamFile(reads=[read])
    seq_parts = gather_reads_for_single_variant(
        samfile=samfile,
        chromosome=chromosome,
        variant=variant)
    print(seq_parts)
    assert len(seq_parts) == 1
    eq_(seq_parts[0].prefix, "ACC")
    eq_(seq_parts[0].variant, "G")
    eq_(seq_parts[0].suffix, "TG")

def test_partitioned_read_sequences_insertion():
    """
    test_partitioned_read_sequences_insertion : Test that read gets correctly
    partitioned for chr1:4 T>TG
    where the sequence for chr1 is assumed to be "ACCTTG"
    """
    # chr1_seq = "ACCTTG"
    chromosome = "chromosome"
    location = 4
    ref = "T"
    alt = "TG"
    variant = Variant(
        chromosome, location, ref, alt, normalize_contig_name=False)

    read = make_read(seq="ACCTGTG", cigar="4M1I2M", mdtag="6")

    samfile = DummySamFile(reads=[read])
    seq_parts = gather_reads_for_single_variant(
        samfile=samfile,
        chromosome=chromosome,
        variant=variant)
    print(seq_parts)
    assert len(seq_parts) == 1
    eq_(seq_parts[0].prefix, "ACCT")
    eq_(seq_parts[0].variant, "G")
    eq_(seq_parts[0].suffix, "TG")

def test_partitioned_read_sequences_deletion():
    """
    test_partitioned_read_sequences_deletion : Test that read gets correctly
    partitioned for chr1:4 TT>T where the sequence for chr1 is assumed to
    be "ACCTTG"
    """
    # chr1_seq = "ACCTTG"
    chromosome = "chromosome"
    location = 4
    ref = "TT"
    alt = "T"
    variant = Variant(
        chromosome, location, ref, alt, normalize_contig_name=False)

    read = make_read(seq="ACCTG", cigar="4M1D1M", mdtag="4^T1")

    samfile = DummySamFile(reads=[read])
    seq_parts = gather_reads_for_single_variant(
        samfile=samfile,
        chromosome=chromosome,
        variant=variant)
    print(seq_parts)
    assert len(seq_parts) == 1
    eq_(seq_parts[0].prefix, "ACCT")
    eq_(seq_parts[0].variant, "")
    eq_(seq_parts[0].suffix, "G")
