"""CIGAR operations, not missing aligned-pair coordinates, define indels."""

import pysam
import pytest

from isovar.allele_read import AlleleRead
from isovar.read_collector import ReadCollector

from .mock_objects import make_pysam_read


@pytest.mark.parametrize("reverse", [False, True])
@pytest.mark.parametrize("md", [None, "8"])
@pytest.mark.parametrize("start,end,accepted", [
    (3, 4, True), (104, 105, True),  # last/first exonic SNPs
    (4, 5, False), (103, 104, False), (4, 104, False),
    (3, 105, False), (3, 104, False), (4, 105, False),
    (3, 3, True), (105, 105, True),  # exonic insertion boundaries
    (4, 4, False), (50, 50, False), (104, 104, False),
])
def test_skip_boundaries_in_direct_and_indexed_collection(tmp_path, reverse, md, start, end, accepted):
    read = make_pysam_read("ACGTACGT", "4M100N4M", mdtag=md)
    read.flag = 16 if reverse else 0
    collector = ReadCollector()
    result = collector.locus_read_from_pysam_aligned_segment(read, start, end)
    assert (result is not None) == accepted
    path = tmp_path / "splice.bam"
    with pysam.AlignmentFile(path, "wb", header={"SQ": [{"SN": "chr1", "LN": 1000}]}) as bam:
        bam.write(read)
    pysam.index(str(path))
    with pysam.AlignmentFile(path) as bam:
        assert bool(collector.get_locus_reads(bam, "chr1", start, end)) == accepted


@pytest.mark.parametrize("reverse", [False, True])
@pytest.mark.parametrize("cigar,seq,md,position,ref,alt,expected", [
    ("3M100N3M1D3M", "CCCAAACCC", "6^A3", 104, "A", "", (3, 3)),
    ("3M1D3M100N3M", "AAACCCGGG", "3^A6", 1, "A", "", (0, 0)),
    ("3M100N3M1I3M", "CCCAAAACCC", "9", 103, "", "A", (3, 4)),
    ("3M1I3M100N3M", "AAAACCCGGG", "9", 1, "", "A", (1, 2)),
    ("3M100N1D3M", "CCCAAA", "3^A3", 104, "A", "", (3, 3)),
    ("3M1D100N3M", "CCCAAA", "3^A3", 4, "A", "", (3, 3)),
    ("3M1D100N1D3M", "CCCAAA", "3^A0^A3", 105, "A", "", (3, 3)),
    ("3M1D1D3M", "AAACCC", "3^AA3", 1, "AA", "", (0, 0)),
    ("3M1I1I3M", "AAAAACCC", "6", 1, "", "AA", (1, 3)),
    ("2=1X1D3=", "ACGCCC", "2A0^A3", 3, "A", "", (2, 2)),
])
def test_normalization_keeps_real_events_and_query_quality_indices(reverse, cigar, seq, md, position, ref, alt, expected):
    read = make_pysam_read(seq, cigar, mdtag=md)
    read.flag = 16 if reverse else 0
    read.query_qualities = list(range(10, 10 + len(seq)))
    assert ReadCollector._left_aligned_indel_interval_for_variant(read, position, ref, alt) == expected
    # The raw read is not rewritten when the canonical query split changes.
    assert list(read.query_qualities) == list(range(10, 10 + len(seq)))


@pytest.mark.parametrize("cigar,seq,md,position,ref,alt", [
    ("3M100N3M1D3M", "AAAAAAAAA", "6^A3", 103, "A", ""),
    ("3M100N3M1I3M", "AAAAAAAAAA", "9", 102, "", "A"),
    ("3M1I2M1D3M", "AAAAAACCC", "5^A3", 3, "A", ""),
    ("3M1D2M1I3M", "AAAAAACCC", "3^A5", 3, "", "A"),
    ("3M2S", "AAAAA", "3", 3, "", "AA"),  # S is not I
    ("2S3M", "AAAAA", "3", 1, "", "AA"),
])
def test_normalization_cannot_jump_discontinuities_or_invent_insertions(cigar, seq, md, position, ref, alt):
    read = make_pysam_read(seq, cigar, mdtag=md)
    assert ReadCollector._left_aligned_indel_interval_for_variant(read, position, ref, alt) is None


@pytest.mark.parametrize("md", [None, "3^A3"])
def test_genuine_deletion_still_extracts_with_or_without_md(md):
    read = make_pysam_read("CCCGGG", "3M1D3M", mdtag=md)
    locus = ReadCollector().locus_read_from_pysam_aligned_segment(
        read, 3, 4, trimmed_base1_start=4, trimmed_ref="A", trimmed_alt="")
    assert AlleleRead.from_locus_read(locus).allele == ""


def test_normalization_does_not_expand_long_introns():
    class NoExpandedPairs:
        def __init__(self, read):
            self.read = read

        def __getattr__(self, name):
            return getattr(self.read, name)

        def get_aligned_pairs(self, *args, **kwargs):
            pytest.fail("Normalization expanded the skipped intron")

    read = NoExpandedPairs(make_pysam_read("CCCAAACCC", "3M10000000N3M1D3M", mdtag="6^A3"))
    assert ReadCollector._left_aligned_indel_interval_for_variant(read, 10000004, "A", "") == (3, 3)


def test_normalization_rejects_inconsistent_md_reference_length():
    class InconsistentReference:
        def __init__(self, read):
            self.read = read

        def __getattr__(self, name):
            return getattr(self.read, name)

        def get_reference_sequence(self):
            return "A"  # not the 7 reference bases consumed by M/D

    read = InconsistentReference(make_pysam_read("AAACCC", "3M1D3M", mdtag="3^A3"))
    assert ReadCollector._left_aligned_indel_interval_for_variant(read, 1, "A", "") is None


def test_normalization_without_query_sequence_is_not_an_exception():
    read = make_pysam_read("AAACCC", "3M1D3M", mdtag="3^A3")
    read.query_sequence = None
    assert ReadCollector._left_aligned_indel_interval_for_variant(read, 1, "A", "") is None
