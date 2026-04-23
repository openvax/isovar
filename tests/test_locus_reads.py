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

from unittest.mock import MagicMock

from varcode import Variant

from isovar.locus_read import LocusRead
from isovar.dataframe_helpers import locus_reads_dataframe
from isovar.read_collector import ReadCollector

from .mock_objects import MockAlignmentFile, make_pysam_read
from .testing_helpers import assert_equal_fields, load_bam, data_path
from .common import eq_
from .genomes_for_testing import grch38


class TrackingReferencePositions(list):
    """
    Count linear index scans over reference positions.
    """

    def __init__(self, values):
        super().__init__(values)
        self.index_calls = 0

    def index(self, *args, **kwargs):
        self.index_calls += 1
        return super().index(*args, **kwargs)


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
    # We know from inspecting the file that one fetched read is spliced across
    # this locus and therefore does not actually overlap the queried interval.
    # The widened fetch window pulls it in, and get_overlap() correctly drops it.
    n_reads_expected -= 1

    print("Found %d sequences in %s" % (n_reads_expected, sam_path_single_variant))
    read_collector = ReadCollector(merge_overlapping_fragments=True)
    reads = read_collector.get_locus_reads(
        sam_all_variants,
        "chr4",
        45802538,
        45802539,
    )
    eq_(sum(read.source_read_count for read in reads), n_reads_expected)

    df = locus_reads_dataframe(
        alignments=sam_all_variants,
        chromosome="chr4",
        base0_start=45802538,
        base0_end=45802539)
    print(df)
    eq_(len(df), n_reads_expected)


def test_get_locus_reads_merges_overlapping_paired_reads():
    """
    Overlapping paired-end reads from the same fragment should collapse into a
    single longer LocusRead while retaining the raw source read count.

    Regression test for GitHub issue #59.
    """
    left_mate = make_pysam_read(
        seq="ACCGTG",
        cigar="6M",
        name="fragment-1",
        reference_start=0,
    )
    right_mate = make_pysam_read(
        seq="CGTGAA",
        cigar="6M",
        name="fragment-1",
        reference_start=2,
    )
    read_collector = ReadCollector(merge_overlapping_fragments=True)
    reads = read_collector.get_locus_reads(
        MockAlignmentFile(references=("chromosome",), reads=[left_mate, right_mate]),
        "chromosome",
        base0_start_inclusive=3,
        base0_end_exclusive=4,
    )
    eq_(len(reads), 1)
    merged_read = reads[0]
    eq_(merged_read.sequence, "ACCGTGAA")
    eq_(merged_read.reference_positions, list(range(8)))
    eq_(merged_read.read_base0_start_inclusive, 3)
    eq_(merged_read.read_base0_end_exclusive, 4)
    eq_(merged_read.source_read_count, 2)


def test_locus_reads_clamps_fetch_start_at_contig_boundary():
    """
    Querying a locus at the first aligned base of a contig should not widen the
    fetch interval to a negative start coordinate.

    Regression test for GitHub issue #161.
    """
    samfile = load_bam("data/primary.chr1.bam")
    read_creator = ReadCollector()
    reads = read_creator.get_locus_reads(
        samfile,
        "chr1",
        base0_start_inclusive=0,
        base0_end_exclusive=1,
    )
    assert isinstance(reads, list)


def test_locus_read_avoids_linear_index_scans_over_reference_positions():
    """
    Regression test for GitHub issue #163.

    Converting the read/reference overlap interval should use a single
    reference-position lookup table instead of repeated linear `list.index`
    scans over the read's reference positions.
    """
    real_read = make_pysam_read(seq="ACCGTG", cigar="6M", mdtag="3G2")
    tracking_positions = TrackingReferencePositions([0, 1, 2, 3, 4, 5])

    mock_read = MagicMock()
    mock_read.query_name = real_read.query_name
    mock_read.query_sequence = real_read.query_sequence
    mock_read.query_qualities = real_read.query_qualities
    mock_read.query_alignment_start = real_read.query_alignment_start
    mock_read.query_alignment_end = real_read.query_alignment_end
    mock_read.is_secondary = False
    mock_read.is_duplicate = False
    mock_read.is_unmapped = False
    mock_read.mapping_quality = real_read.mapping_quality
    mock_read.get_reference_positions.return_value = tracking_positions

    read_collector = ReadCollector()
    result = read_collector.locus_read_from_pysam_aligned_segment(
        mock_read,
        base0_start_inclusive=3,
        base0_end_exclusive=4,
    )
    assert result is not None
    eq_(tracking_positions.index_calls, 0)


def _make_mock_pysam_read_with_none_mapq():
    """
    Create a mock pysam read whose mapping_quality is None.

    pysam.AlignedSegment does not allow setting mapping_quality to None
    directly, so we use a MagicMock instead.
    """
    real_read = make_pysam_read(seq="ACCGTG", cigar="6M", mdtag="3G2")
    mock_read = MagicMock()
    mock_read.query_name = real_read.query_name
    mock_read.query_sequence = real_read.query_sequence
    mock_read.query_qualities = real_read.query_qualities
    mock_read.is_secondary = False
    mock_read.is_duplicate = False
    mock_read.is_unmapped = False
    mock_read.get_reference_positions.return_value = list(range(6))
    mock_read.mapping_quality = None
    return mock_read


def test_locus_read_none_mapq_with_min_mapping_quality_zero():
    """
    When min_mapping_quality=0, a read with mapping_quality=None should NOT
    be skipped (0 means accept everything). The None should be treated as 0.

    Regression test for GitHub issue #127.
    """
    mock_read = _make_mock_pysam_read_with_none_mapq()
    read_collector = ReadCollector(min_mapping_quality=0)
    result = read_collector.locus_read_from_pysam_aligned_segment(
        mock_read,
        base0_start_inclusive=3,
        base0_end_exclusive=4,
    )
    assert result is not None, (
        "Expected read with mapping_quality=None to be accepted when "
        "min_mapping_quality=0, but it was skipped"
    )


def test_locus_read_none_mapq_with_min_mapping_quality_nonzero():
    """
    When min_mapping_quality > 0, a read with mapping_quality=None should
    be skipped.

    Regression test for GitHub issue #127.
    """
    mock_read = _make_mock_pysam_read_with_none_mapq()
    read_collector = ReadCollector(min_mapping_quality=1)
    result = read_collector.locus_read_from_pysam_aligned_segment(
        mock_read,
        base0_start_inclusive=3,
        base0_end_exclusive=4,
    )
    assert result is None, (
        "Expected read with mapping_quality=None to be skipped when "
        "min_mapping_quality=1, but it was accepted"
    )


def test_locus_read_insertion_locus_no_insertion_in_read():
    """
    When an insertion variant is queried (reference_interval_size == 0) but the
    read does NOT carry the insertion, the read positions flanking the locus
    should be adjacent (differ by 1) and the returned LocusRead should have
    an empty read interval (start == end).

    Regression test for GitHub issue #126 — previously the code compared
    read_base0_after_insertion to itself instead of read_base0_before_insertion,
    making the adjacent-positions branch dead.
    """
    # Reference: ACCTTG (positions 0-5). Insertion locus between pos 3 and 4.
    # This read has NO insertion — reference positions are contiguous.
    pysam_read = make_pysam_read(seq="ACCTTG", cigar="6M", mdtag="6")
    read_collector = ReadCollector()
    result = read_collector.locus_read_from_pysam_aligned_segment(
        pysam_read,
        base0_start_inclusive=4,
        base0_end_exclusive=4,
    )
    assert result is not None, "Read without insertion should still be returned"
    eq_(result.read_base0_start_inclusive, 4)
    eq_(result.read_base0_end_exclusive, 4)


def test_locus_read_insertion_locus_with_insertion_in_read():
    """
    When an insertion variant is queried (reference_interval_size == 0) and the
    read DOES carry the insertion, the returned LocusRead should have a
    non-empty read interval spanning the inserted bases.

    Companion to the test above to verify the else-branch still works.
    """
    # Reference: ACCTTG. Insertion of G between pos 3 and 4.
    # Read: ACCTGTG — 4M1I2M — positions [0,1,2,3,None,4,5]
    pysam_read = make_pysam_read(seq="ACCTGTG", cigar="4M1I2M", mdtag="6")
    read_collector = ReadCollector()
    result = read_collector.locus_read_from_pysam_aligned_segment(
        pysam_read,
        base0_start_inclusive=4,
        base0_end_exclusive=4,
    )
    assert result is not None, "Read with insertion should be returned"
    eq_(result.read_base0_start_inclusive, 4)
    eq_(result.read_base0_end_exclusive, 5)


def test_locus_read_insertion_locus_with_insertion_at_start_of_read():
    """
    An insertion-supporting read can begin at the insertion locus and still
    provide a non-empty inserted interval.

    Regression test for GitHub issue #49.
    """
    pysam_read = make_pysam_read(
        seq="GT",
        cigar="1I1M",
        mdtag="1",
        reference_start=3,
    )
    read_collector = ReadCollector()
    result = read_collector.locus_read_from_pysam_aligned_segment(
        pysam_read,
        base0_start_inclusive=3,
        base0_end_exclusive=3,
    )
    assert result is not None, "Insertion at the start of a read should be kept"
    eq_(result.read_base0_start_inclusive, 0)
    eq_(result.read_base0_end_exclusive, 1)


def test_locus_read_insertion_locus_with_insertion_at_end_of_read():
    """
    An insertion-supporting read can end at the insertion locus and still
    provide a non-empty inserted interval.

    Regression test for GitHub issue #49.
    """
    pysam_read = make_pysam_read(
        seq="AG",
        cigar="1M1I",
        mdtag="1",
        reference_start=2,
    )
    read_collector = ReadCollector()
    result = read_collector.locus_read_from_pysam_aligned_segment(
        pysam_read,
        base0_start_inclusive=3,
        base0_end_exclusive=3,
    )
    assert result is not None, "Insertion at the end of a read should be kept"
    eq_(result.read_base0_start_inclusive, 1)
    eq_(result.read_base0_end_exclusive, 2)


def test_locus_read_left_aligns_equivalent_homopolymer_insertion():
    """
    A read whose insertion is aligned to the right within a homopolymer should
    still be rewritten to the canonical left-aligned read interval for the
    queried variant.

    Regression test for GitHub issue #79.
    """
    variant = Variant(
        "1",
        1,
        "A",
        "AA",
        grch38,
        normalize_contig_names=False,
    )
    pysam_read = make_pysam_read(
        seq="AAAA",
        cigar="2M1I1M",
        mdtag="3",
        reference_start=0,
    )
    read_collector = ReadCollector()
    reads = read_collector.locus_reads_overlapping_variant(
        MockAlignmentFile(references=("1",), reads=[pysam_read]),
        variant,
        chromosome="1",
    )
    assert len(reads) == 1
    eq_(reads[0].read_base0_start_inclusive, 1)
    eq_(reads[0].read_base0_end_exclusive, 2)


def test_locus_read_left_aligns_equivalent_homopolymer_deletion():
    """
    A deletion aligned to the right within a homopolymer should also be
    rewritten to the canonical left-aligned query interval.
    """
    variant = Variant(
        "1",
        1,
        "AA",
        "A",
        grch38,
        normalize_contig_names=False,
    )
    pysam_read = make_pysam_read(
        seq="AAA",
        cigar="2M1D1M",
        mdtag="2^A1",
        reference_start=0,
    )
    read_collector = ReadCollector()
    reads = read_collector.locus_reads_overlapping_variant(
        MockAlignmentFile(references=("1",), reads=[pysam_read]),
        variant,
        chromosome="1",
    )
    assert len(reads) == 1
    eq_(reads[0].read_base0_start_inclusive, 1)
    eq_(reads[0].read_base0_end_exclusive, 1)
