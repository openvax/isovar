"""Audit the recorded full-depth #234 coordinate-sharing evidence without huge BAMs."""

import gzip
import json
from pathlib import Path

import pytest

from tests.data.osteosarc.expansion.inventory import digest
from tests.evidence_audit_helpers import (
    COLLECTION_PAIR_IDENTITY_KEYS, COLLECTION_RESULT_KEYS, PIPELINE_PAIR_IDENTITY_KEYS, PIPELINE_RESULT_KEYS,
    changed_sources, load_pinned_package)


EXPANSION = Path(__file__).parent / "data/osteosarc/expansion"
RESULTS = EXPANSION / "performance/coordinate-sharing-results"
COLLECTION_LABELS = {"t1", "t1primary", "t3dedup", "t3dedupprimary", "t3tagged", "ontdync", "ontdel", "bulk",
                     "t3exoc4", "profile_t1"}
# Recorded digests (COORDINATE_SHARING.md, "Evidence"): they pin what produced
# this package without freezing the shipped scripts, which may keep evolving.
COLLECTION_BENCHMARK_SHA256 = "6429ffaa26a01d3f28f6d5324ae123f33ac85ec26a050a49be414fb84c4c883e"
PIPELINE_BENCHMARK_SHA256 = "1492afb2987b28f01ff3c864b4b97c9e7bde514cd93a3a91a1d0133ef9c3792e"
SCANNER_SHA256 = "d08d08732646ea905aea5e5f3622d8382a2c7ae477cdd9fbbec0bc053308afac"
BOUND_SCAN_SHA256 = "c9f789b24980f77c5cfd6bd74fdc26032408e5201dac02d3bbf49119d6cc2457"


@pytest.fixture(scope="module")
def evidence():
    return load_pinned_package(RESULTS, EXPANSION / "collection_report.py")


@pytest.fixture(scope="module")
def bound_scan():
    path = RESULTS / "bound-scan.json.gz"
    assert digest(path) == BOUND_SCAN_SHA256
    return json.loads(gzip.decompress(path.read_bytes()))


def test_every_run_records_the_pinned_measuring_scripts(evidence):
    for run in evidence["collection"]["runs"].values():
        assert run["identity"]["benchmark_sha256"] == COLLECTION_BENCHMARK_SHA256
    for run in evidence["pipeline"]["runs"].values():
        assert run["identity"]["benchmark_sha256"] == PIPELINE_BENCHMARK_SHA256


def test_per_call_sharing_preserves_every_full_depth_collection(evidence):
    data = evidence["collection"]
    assert len(data["comparisons"]) == 10
    assert len({c["after"] for c in data["comparisons"]}) == 10
    assert sum(c["allocation_instrumented"] for c in data["comparisons"]) == 1
    assert {c["after"].removeprefix("final_") for c in data["comparisons"]} == COLLECTION_LABELS
    modes = set()
    for comparison in data["comparisons"]:
        assert comparison["exact_locus_and_allele_reads"] is True
        before, after = (data["runs"][comparison[key]] for key in ("before", "after"))
        assert (before["identity"]["isovar"], after["identity"]["isovar"]) == ("1.8.3", "1.8.4")
        assert before["result"]["status"] == after["result"]["status"] == "ok"
        for key in COLLECTION_PAIR_IDENTITY_KEYS:
            assert before["identity"][key] == after["identity"][key]
        assert changed_sources(before, after) == {"isovar/__init__.py", "isovar/read_collector.py"}
        for key in COLLECTION_RESULT_KEYS:
            assert before["result"][key] == after["result"][key]
        modes.add(after["identity"]["mode"])
    assert modes == {"defaults", "primary_only"}


def test_t3_tagged_pipeline_ranked_proteins_and_rna_are_unchanged(evidence):
    data = evidence["pipeline"]
    # Two same-session pairs with the run order reversed: a single whole-process
    # RSS observation on this host is dominated by host memory state.
    assert sorted((c["before"], c["after"]) for c in data["comparisons"]) == [
        ("baseline_t3tagged", "final_t3tagged"), ("baseline_t3tagged_repeat", "final_t3tagged_repeat")]
    for comparison in data["comparisons"]:
        assert comparison["exact_rna_and_ranked_proteins"] is True
        before, after = (data["runs"][comparison[key]] for key in ("before", "after"))
        assert (before["identity"]["isovar"], after["identity"]["isovar"]) == ("1.8.3", "1.8.4")
        for key in PIPELINE_PAIR_IDENTITY_KEYS:
            assert before["identity"][key] == after["identity"][key]
        assert changed_sources(before, after) == {
            "isovar/__init__.py", "isovar/read_collector.py", "tests/data/osteosarc/expansion/collection_benchmark.py"}
        for run in (before, after):
            assert run["result"]["status"] == run["result"]["validation_status"] == "ok"
            assert all(m["status"] == "ok" for m in run["measurements"].values())
        for key in PIPELINE_RESULT_KEYS:
            assert before["result"][key] == after["result"][key]
        assert after["result"]["checked_protein_count"] == 4351


def test_frozen_collection_observations_keep_memory_and_show_no_cpu_regression(evidence):
    # These check the recorded observations, not machine-dependent thresholds
    # for code executed during CI. No speedup is claimed.
    runs = evidence["collection"]["runs"]
    totals = dict(baseline=0.0, final=0.0)
    for label in COLLECTION_LABELS - {"profile_t1"}:
        before, after = (runs[f"{phase}_{label}"]["result"] for phase in ("baseline", "final"))
        assert abs(after["locus_collection_peak_rss_bytes"] / before["locus_collection_peak_rss_bytes"] - 1) < 0.03
        totals["baseline"] += before["locus_collection_cpu_seconds"]
        totals["final"] += after["locus_collection_cpu_seconds"]
    assert totals["final"] <= totals["baseline"]


def test_allocation_profile_shows_no_added_retained_memory(evidence):
    runs = evidence["collection"]["runs"]
    before, after = (runs[f"{phase}_profile_t1"]["allocations"] for phase in ("baseline", "final"))
    assert before["stage"] == after["stage"] == "before_mate_merging"
    assert before["locus_read_count"] == after["locus_read_count"] == 34071
    assert before["query_bases"] == after["query_bases"] == 30619733
    assert after["traced_current_bytes"] <= before["traced_current_bytes"]
    for profile in (before, after):
        assert profile["traced_peak_bytes"] >= profile["traced_current_bytes"] > 0
        assert profile["tracer_overhead_bytes"] > 0
        assert profile["locations"]


def test_bound_scan_is_tied_to_the_packaged_collection_evidence(bound_scan, evidence):
    identity, reproduction = bound_scan["identity"], bound_scan["reproduction"]
    assert identity["isovar"] == "1.8.3"
    assert identity["lru_maxsize"] == [65536]
    assert identity["scanner_sha256"] == reproduction["scanner_sha256"] == SCANNER_SHA256
    exoc4 = evidence["collection"]["runs"]["final_t3exoc4"]
    assert identity["inventory_sha256"] == exoc4["identity"]["inventory_sha256"]
    assert identity["receipt_bam_sha256"][exoc4["identity"]["source_id"]] == exoc4["identity"]["source_bam_sha256"]
    scanned = {row["source_id"] for row in bound_scan["rows"]}
    assert len(scanned) == 134 and set(identity["receipt_bam_sha256"]) == scanned
    # Some scanned sources' full BAMs later left the local cache; every row whose
    # source remained was reproduced exactly by the committed scanner.
    assert reproduction["identical_rows"] == reproduction["rows"] == 3009
    assert reproduction["sources"] == 88
    assert len(reproduction["sources_without_bam"]) == 46 and set(reproduction["sources_without_bam"]) <= scanned


def test_bound_scan_supports_the_documented_table_size_claims(bound_scan, evidence):
    rows, int_bytes = bound_scan["rows"], bound_scan["identity"]["int_bytes"]
    assert len(rows) == 4563
    assert max(row["distinct"] for row in rows) < 65536
    largest = max(rows, key=lambda row: row["distinct"])
    assert (largest["source_id"], largest["variant_id"]) == ("58c4d68e5f7e55b5", "EXOC4-chr7-133274996")
    assert (largest["reads"], largest["lookups"], largest["distinct"]) == (2370, 3326558, 46702)
    assert largest["reads"] == evidence["collection"]["runs"]["final_t3exoc4"]["result"]["locus_read_count"]
    assert round(largest["hits"] / largest["lookups"], 3) == 0.986
    assert largest["table_bytes"] == 2621528
    assert (largest["lookups"] - largest["distinct"]) * int_bytes == 91835968
    large = [row["table_bytes"] / row["distinct"] for row in rows if row["distinct"] >= 10000]
    assert len(large) == 47
    assert (round(min(large), 1), round(max(large), 1)) == (27.1, 56.4)
    unpaid = [row for row in rows if row["table_bytes"] > (row["lookups"] - row["distinct"]) * int_bytes]
    assert len(unpaid) == 457
    worst = max(unpaid, key=lambda row: row["table_bytes"])
    assert (worst["source_id"], worst["variant_id"], worst["reads"], worst["table_bytes"]) == (
        "53f498a544883d51", "WWC1-chr5-168399513", 7, 294992)
