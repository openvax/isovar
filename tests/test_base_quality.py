"""Reported QUAL survives extraction; optional assessments never rewrite support."""

from copy import copy
import pickle

import pytest

from isovar.allele_read import AlleleRead
from isovar.base_quality import BaseQualityPolicy
from isovar.read_collector import ReadCollector, _CompactLocusRead
from isovar.rna_evidence import EvidenceSets
from .mock_objects import make_pysam_read


def observation(sequence="ACGT", cigar="4M", start=1, end=2, qualities=(30, 20, 10, 0),
                flag=0, **kwargs):
    read = make_pysam_read(sequence, cigar, name="original", mdtag=kwargs.pop("md", None))
    read.flag = flag
    read.query_qualities = qualities
    locus = ReadCollector(**kwargs).locus_read_from_pysam_aligned_segment(read, start, end)
    return locus, AlleleRead.from_locus_read(locus)


@pytest.mark.parametrize("threshold", [-1, 255, 1.5, True, None, "20"])
def test_invalid_thresholds(threshold):
    with pytest.raises(ValueError, match="integer"):
        BaseQualityPolicy(threshold)


@pytest.mark.parametrize("qualities,threshold,status", [
    ((20, 20, 20, 20), 20, "passed"), ((20, 20, 20, 20), 21, "failed"),
    ((0, 0, 0, 0), 0, "passed"), ((0, 0, 0, 0), 1, "failed"),
    (None, 0, "unassessed"), (None, 20, "unassessed"),
    ((40, 255, 40, 40), 20, "unassessed"),
])
def test_constant_missing_zero_and_boundary_qualities(qualities, threshold, status):
    _, read = observation(qualities=qualities)
    assessment, = BaseQualityPolicy(threshold).assess(read)
    assert assessment.status == status
    assert read.allele == "C"
    assert BaseQualityPolicy(threshold).description()["calibrated_probability"] is False


@pytest.mark.parametrize("flag", [0, 16])
@pytest.mark.parametrize("seq,cigar,start,end,indices,role,allele", [
    ("ACGT", "4M", 1, 3, (1, 2), "allele_bases", "CG"),
    ("ACGT", "1M2I1M", 1, 1, (1, 2), "allele_bases", "CG"),
    ("ACGT", "2M2D2M", 2, 4, (1, 2), "flanking_anchors", ""),
    ("ACGT", "4M", 2, 2, (1, 2), "flanking_anchors", ""),
])
def test_original_query_indices_and_compact_path(flag, seq, cigar, start, end, indices, role, allele):
    locus, read = observation(seq, cigar, start, end, flag=flag)
    compact = AlleleRead.from_locus_read(_CompactLocusRead.from_locus_read(locus))
    for observed in (read, compact, read.with_compatible_transcript_ids(["tx"])):
        assert observed.allele == allele
        footprint, = observed.source_allele_qualities
        assert footprint.quality_indices == indices
        assert footprint.quality_scores == (20, 10)
        assert footprint.role == role
        assert footprint.source_alignment[1][-1] == bool(flag & 16)
        assert observed.quality_scores == (30, 20, 10, 0)
        assert BaseQualityPolicy(10).assess(observed)[0].status == "passed"
        assert BaseQualityPolicy(11).assess(observed)[0].status == "failed"


def test_trimming_preserves_full_vector_and_original_coordinates():
    locus, read = observation("NTACGTAN", "2S4M2S", 1, 2, qualities=tuple(range(8)))
    assert read.sequence == "ACGT" and read.quality_scores == (2, 3, 4, 5)
    footprint, = read.source_allele_qualities
    assert footprint.query_interval == footprint.canonical_query_interval == (3, 4)
    assert footprint.quality_indices == (3,) and footprint.quality_scores == (3,)
    _, with_clips = observation("NTACGTAN", "2S4M2S", 1, 2, qualities=tuple(range(8)),
                               use_soft_clipped_bases=True)
    assert with_clips.sequence == "TACGTA" and with_clips.quality_scores == (1, 2, 3, 4, 5, 6)
    assert with_clips.source_allele_qualities == read.source_allele_qualities


@pytest.mark.parametrize("flag", [0, 16])
def test_repeat_shifted_deletion_uses_original_anchors(flag):
    record = make_pysam_read("CCCAAACCC", "6M1D3M", mdtag="6^A3")
    record.flag = flag
    record.query_qualities = [10, 10, 10, 10, 10, 30, 35, 10, 10]
    locus = ReadCollector().locus_read_from_pysam_aligned_segment(
        record, 3, 4, trimmed_base1_start=4, trimmed_ref="A", trimmed_alt="")
    read = AlleleRead.from_locus_read(locus)
    footprint, = read.source_allele_qualities
    assert footprint.canonical_query_interval == (3, 3)
    assert footprint.query_interval == (6, 6)
    assert footprint.quality_indices == (5, 6) and footprint.quality_scores == (30, 35)
    assert BaseQualityPolicy(30).assess(read)[0].status == "passed"


def test_terminal_deletion_does_not_borrow_a_soft_clip_quality():
    record = make_pysam_read("AACCC", "2S1D3M", mdtag="0^A3")
    record.query_qualities = [40] * 5
    locus = ReadCollector().locus_read_from_pysam_aligned_segment(
        record, 0, 1, trimmed_base1_start=1, trimmed_ref="A", trimmed_alt="")
    read = AlleleRead.from_locus_read(locus)
    footprint, = read.source_allele_qualities
    assert footprint.quality_indices == (None, 2)
    assert footprint.quality_scores == (None, 40)
    assert BaseQualityPolicy(20).assess(read)[0].status == "unassessed"
    read.source_allele_qualities = (footprint._replace(quality_scores=(None, 10)),)
    assert BaseQualityPolicy(20).assess(read)[0].status == "failed"


@pytest.mark.parametrize("second_allele,second_quality,expected", [("C", 10, "failed"), ("A", 10, "unassessed"),
                                                                  ("C", None, "unassessed"), ("C", 40, "passed")])
def test_mate_qualities_are_assessed_individually(second_allele, second_quality, expected):
    first, _ = observation(flag=1 | 64, qualities=[30] * 4)
    second, _ = observation(sequence="A" + second_allele + "GT", flag=1 | 128,
                            qualities=None if second_quality is None else [second_quality] * 4)
    merged = ReadCollector._merge_locus_read_pair(first, second)
    assert merged is not None
    read = AlleleRead.from_locus_read(merged)
    assert read.allele == "C" and read.source_read_count == 2
    assessments = BaseQualityPolicy(20).assess(read)
    assert [a.status for a in assessments] == ["passed", expected]
    assert len(read.source_allele_qualities) == 2
    support = EvidenceSets(["sample", "input"], BaseQualityPolicy(20)).support([read])
    assert (support["reads"], support["fragments"]) == (2, 1)
    assert support["base_quality"]["passed"]["reads"] == (2 if expected == "passed" else 1)
    assert sum(support["base_quality"][s]["reads"] for s in ("passed", "failed", "unassessed")) == 2


def test_alternative_placements_do_not_count_or_qualify_twice():
    _, high = observation()
    low = copy(high)
    original, = high.source_allele_qualities
    key, placement = original.source_alignment
    alternative = (key, (placement[0], placement[1] + 100, *placement[2:]))
    low.source_alignments = (alternative,)
    low.source_allele_qualities = (original._replace(source_alignment=alternative, quality_scores=(0,)),)
    support = EvidenceSets(["s", "input"], BaseQualityPolicy(20)).support([high, high, low])
    assert support["reads"] == support["base_quality"]["failed"]["reads"] == 1
    assert support["base_quality"]["passed"]["reads"] == 0
    low.source_allele_qualities = ()
    support = EvidenceSets(["s", "input"], BaseQualityPolicy(20)).support([high, low])
    assert support["base_quality"]["unassessed"]["reads"] == 1


def test_legacy_identity_and_metadata_roundtrip():
    _, read = observation()
    missing = copy(read)
    del missing.quality_scores
    del missing.source_allele_qualities
    assert read == missing and hash(read) == hash(missing)
    assert BaseQualityPolicy(20).assess(missing)[0].reason == "missing_provenance"
    cloned = missing.with_compatible_transcript_ids(["tx"])
    assert cloned.quality_scores is None and cloned.source_allele_qualities == ()
    restored = pickle.loads(pickle.dumps(read))
    assert restored.quality_scores == read.quality_scores
    assert restored.source_allele_qualities == read.source_allele_qualities
    legacy = AlleleRead("A", "C", "GT", "legacy", source_read_count=2)
    support = EvidenceSets(["s", "input"], BaseQualityPolicy(20)).support([legacy])
    assert support["reads"] == support["base_quality"]["unassessed"]["reads"] == 2
    assert support["base_quality"]["unassessed"]["evidence_set_id"] is None
    with pytest.raises(ValueError, match="length"):
        AlleleRead("A", "C", "GT", "bad", quality_scores=[30])
