# Copyright (c) 2016-2018. Mount Sinai School of Medicine
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
from isovar.variant_reads import reads_supporting_variant
from isovar.allele_reads import AlleleRead

from mock_read_data import DummySamFile, make_read
from nose.tools import eq_


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
    variant_reads = reads_supporting_variant(
        samfile=samfile,
        chromosome=chromosome,
        variant=variant)
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
    chromosome = "chromosome"
    location = 4
    ref = "T"
    alt = "TG"
    variant = Variant(
        chromosome, location, ref, alt, normalize_contig_name=False)

    read = make_read(seq=b"ACCTGTG", cigar="4M1I2M", mdtag="6")

    samfile = DummySamFile(reads=[read])
    variant_reads = reads_supporting_variant(
        samfile=samfile,
        chromosome=chromosome,
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
    chromosome = "chromosome"
    location = 4
    ref = "TT"
    alt = "T"
    variant = Variant(
        chromosome, location, ref, alt, normalize_contig_name=False)

    read = make_read(seq="ACCTG", cigar="4M1D1M", mdtag="4^T1")

    samfile = DummySamFile(reads=[read])
    variant_reads = reads_supporting_variant(
        samfile=samfile,
        chromosome=chromosome,
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
