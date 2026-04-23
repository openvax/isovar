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

from .common import eq_ 
from varcode import Variant

from isovar.allele_read import AlleleRead
from isovar.read_collector import ReadCollector

from .mock_objects import MockAlignmentFile, make_pysam_read
from .genomes_for_testing import grch38


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


def test_read_evidence_counts_snv_reads_at_read_boundaries():
    """
    Reads whose SNV locus is the first or last aligned base should still be
    counted.

    Regression test for GitHub issue #55.
    """
    chromosome = "1"
    variant = Variant(
        chromosome,
        4,
        "T",
        "G",
        grch38,
        normalize_contig_names=False)

    reads = [
        make_pysam_read(
            seq="GAAA",
            cigar="4M",
            mdtag="0T3",
            name="alt-start",
            reference_start=3),
        make_pysam_read(
            seq="TAAA",
            cigar="4M",
            mdtag="4",
            name="ref-start",
            reference_start=3),
        make_pysam_read(
            seq="AAAG",
            cigar="4M",
            mdtag="3T0",
            name="alt-end",
            reference_start=0),
        make_pysam_read(
            seq="AAAT",
            cigar="4M",
            mdtag="4",
            name="ref-end",
            reference_start=0),
    ]

    read_creator = ReadCollector()
    read_evidence = read_creator.read_evidence_for_variant(
        variant=variant,
        alignment_file=MockAlignmentFile(references=(chromosome,), reads=reads))

    eq_(read_evidence.alt_read_names, {"alt-start", "alt-end"})
    eq_(read_evidence.ref_read_names, {"ref-start", "ref-end"})


def test_read_evidence_counts_insertion_reads_at_read_boundaries():
    """
    Insertion-supporting reads that begin or end at the insertion locus should
    be counted as alt reads instead of being dropped.

    Regression test for GitHub issue #49. Also serves as an indel-count
    regression for GitHub issue #23.
    """
    chromosome = "1"
    variant = Variant(
        chromosome,
        3,
        "A",
        "AG",
        grch38,
        normalize_contig_names=False)

    reads = [
        make_pysam_read(
            seq="GT",
            cigar="1I1M",
            mdtag="1",
            name="alt-start",
            reference_start=3),
        make_pysam_read(
            seq="T",
            cigar="1M",
            mdtag="1",
            name="ref-start",
            reference_start=3),
        make_pysam_read(
            seq="AG",
            cigar="1M1I",
            mdtag="1",
            name="alt-end",
            reference_start=2),
        make_pysam_read(
            seq="A",
            cigar="1M",
            mdtag="1",
            name="ref-end",
            reference_start=2),
    ]

    read_creator = ReadCollector()
    read_evidence = read_creator.read_evidence_for_variant(
        variant=variant,
        alignment_file=MockAlignmentFile(references=(chromosome,), reads=reads))

    eq_(read_evidence.alt_read_names, {"alt-start", "alt-end"})
    eq_(read_evidence.ref_read_names, {"ref-start", "ref-end"})


def test_partitioned_read_sequences_snv_at_last_exonic_base_before_splice():
    """
    A variant on the last nucleotide of an exon should still be recovered from
    a spliced read.

    Regression test for GitHub issue #24.
    """
    chromosome = "1"
    variant = Variant(
        chromosome,
        4,
        "T",
        "G",
        grch38,
        normalize_contig_names=False)
    read = make_pysam_read(
        seq="ACCGTGGA",
        cigar="4M100N4M",
        name="splice-before",
        reference_start=0)
    read_creator = ReadCollector()
    variant_reads = read_creator.allele_reads_supporting_variant(
        variant=variant,
        alignment_file=MockAlignmentFile(references=(chromosome,), reads=[read]))
    assert len(variant_reads) == 1
    eq_(variant_reads[0].allele, "G")


def test_partitioned_read_sequences_snv_at_first_exonic_base_after_splice():
    """
    A variant on the first nucleotide of an exon should still be recovered from
    a spliced read.

    Regression test for GitHub issue #24.
    """
    chromosome = "1"
    variant = Variant(
        chromosome,
        105,
        "T",
        "G",
        grch38,
        normalize_contig_names=False)
    read = make_pysam_read(
        seq="ACCTGGGA",
        cigar="4M100N4M",
        name="splice-after",
        reference_start=0)
    read_creator = ReadCollector()
    variant_reads = read_creator.allele_reads_supporting_variant(
        variant=variant,
        alignment_file=MockAlignmentFile(references=(chromosome,), reads=[read]))
    assert len(variant_reads) == 1
    eq_(variant_reads[0].allele, "G")


def test_partitioned_read_sequences_left_align_homopolymer_insertion():
    """
    A read whose insertion is aligned to the right within a homopolymer should
    still count as supporting the left-aligned insertion variant.

    Regression test for GitHub issue #79.
    """
    chromosome = "1"
    variant = Variant(
        chromosome,
        1,
        "A",
        "AA",
        grch38,
        normalize_contig_names=False,
    )
    read = make_pysam_read(
        seq="AAAA",
        cigar="2M1I1M",
        mdtag="3",
        reference_start=0,
    )
    read_creator = ReadCollector()
    read_evidence = read_creator.read_evidence_for_variant(
        variant=variant,
        alignment_file=MockAlignmentFile(references=(chromosome,), reads=[read]),
    )
    eq_(read_evidence.alt_read_names, {"dummy"})
    eq_(read_evidence.ref_read_names, set())


def test_partitioned_read_sequences_left_align_homopolymer_deletion():
    """
    A read whose deletion is aligned to the right within a homopolymer should
    still count as supporting the left-aligned deletion variant.
    """
    chromosome = "1"
    variant = Variant(
        chromosome,
        1,
        "AA",
        "A",
        grch38,
        normalize_contig_names=False,
    )
    read = make_pysam_read(
        seq="AAA",
        cigar="2M1D1M",
        mdtag="2^A1",
        reference_start=0,
    )
    read_creator = ReadCollector()
    read_evidence = read_creator.read_evidence_for_variant(
        variant=variant,
        alignment_file=MockAlignmentFile(references=(chromosome,), reads=[read]),
    )
    eq_(read_evidence.alt_read_names, {"dummy"})
    eq_(read_evidence.ref_read_names, set())
