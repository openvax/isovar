"""Independent-product audit: original RNA, not pooled processing copies."""

from collections import Counter

import pysam
import pytest

from tests.data.osteosarc.figure_comparisons import nr2f2, nr2f2_libraries


DATA = nr2f2_libraries.load()
SOURCES = {s["source_id"]: s for s in DATA["sources"]}
CEGAT = "fbeb5d8a415e4d18"


def test_assemblies_and_unassessed_products_are_explicit():
    assert Counter(s["status"] for s in SOURCES.values()) == {
        "ok": 23,
        "not_assessed": 2,
    }
    assert (
        DATA["references"]["GRCh37"]["dna"].upper()
        == DATA["references"]["GRCh38"]["dna"].upper()
    )
    assert DATA["mapping"] == dict(
        chrom="chr15", start=96875510, end=96875600, strand="+", chain_id="16"
    )
    assert DATA["metadata_sha256"] == nr2f2.load()["metadata_sha256"]


@pytest.mark.parametrize(
    "source_id", sorted(k for k, s in SOURCES.items() if s["status"] == "ok")
)
@pytest.mark.parametrize(
    "name,q,m",
    [
        ("q20", 20, 20),
        ("q30", 30, 20),
        ("q10", 10, 20),
        ("mapq1", 20, 1),
        ("sequence_only", None, 20),
    ],
)
def test_original_records_reproduce_each_threshold(source_id, name, q, m):
    source = SOURCES[source_id]
    focal, deletion = (
        (96875527, (96875576, 96875579))
        if source["assembly"] == "GRCh37"
        else (96332298, (96332347, 96332350))
    )
    audit = nr2f2.recount(
        source,
        DATA["references"][source["assembly"]],
        min_quality=q,
        min_mapq=m,
        focal_position=focal,
        deletion_interval=deletion,
    )
    assert audit == source["audits"][name]
    if source_id != CEGAT:
        # Even the explicitly unqualified screen has no exact deletion sequence.
        assert audit["counts"]["deletion"].get("deletion", 0) == 0
        assert not audit["deletion_records"]


def test_known_processing_copies_are_not_independent_evidence():
    def names(key):
        return {line.split("\t")[0] for line in SOURCES[key]["records"]}

    assert names("57acecacfe1b2636") == names("a660879c8d40a735")  # Tempus deliveries
    assert names("c712686597d2160a") == names("c89442609fffd3f1")  # T1 CellRanger
    assert len(names("47b695c1a686c058") & names("261f8be4ce5aa8c8")) == 402  # T1 BG
    assert names("7010dd389c50467d") < names("53f498a544883d51")  # ONT dedup subset


def test_missing_qualities_are_not_silently_promoted_to_q20():
    source = SOURCES["0066232879babe83"]
    header = pysam.AlignmentHeader.from_text(source["header"])
    assert all(
        pysam.AlignedSegment.fromstring(line, header).query_qualities is None
        for line in source["records"]
    )
    assert source["audits"]["q20"]["counts"]["deletion"] == {}
    assert source["audits"]["sequence_only"]["counts"]["deletion"] == {"reference": 19}


def test_cegat_support_stays_one_endpoint_family():
    audit = SOURCES[CEGAT]["audits"]["q20"]
    assert audit["counts"]["deletion"] == {"reference": 23, "deletion": 5}
    assert (
        audit["deletion_alignment_paths"] == audit["deletion_fragment_endpoints"] == 1
    )
