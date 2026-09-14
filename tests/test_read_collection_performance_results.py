"""Audit the recorded full-depth #229 evidence without loading huge BAMs."""

from pathlib import Path

import pytest

from tests.evidence_audit_helpers import (
    COLLECTION_PAIR_IDENTITY_KEYS, COLLECTION_RESULT_KEYS, PIPELINE_RESULT_KEYS, changed_sources, load_pinned_package)


EXPANSION = Path(__file__).parent / "data/osteosarc/expansion"
RESULTS = EXPANSION / "performance/collection-results"
# Recorded measuring-script digests (COLLECTION_MEMORY.md, "Evidence"): the
# ordinary pairs ran collection_benchmark.py from c165e9c, and the instrumented
# pair ran the 1e47bd6 tracer-stop fix, which only changes --allocations runs.
ORDINARY_BENCHMARK_SHA256 = "a3ced87fef699537a88a4c6d908c4131c51c250c647b250bb5b1728e53cdc254"
PROFILE_BENCHMARK_SHA256 = "c2f7cd30df297ff8dad2fb21ede72ee3d5fac23a6d2ff7452706312215f16c87"


@pytest.fixture(scope="module")
def collection_evidence():
    return load_pinned_package(RESULTS, EXPANSION / "collection_report.py")


def test_full_depth_locus_allele_and_evidence_comparisons_are_exact(collection_evidence):
    data = collection_evidence["collection"]
    ordinary = [c for c in data["comparisons"] if not c["allocation_instrumented"]]
    assert len(ordinary) == 8
    assert len(data["comparisons"]) == 9  # Separate instrumented T1 pair.
    modes = set()
    for comparison in data["comparisons"]:
        assert comparison["exact_locus_and_allele_reads"] is True
        before, after = [data["runs"][comparison[k]] for k in ("before", "after")]
        assert before["identity"]["isovar"] == "1.8.2"
        assert after["identity"]["isovar"] == "1.8.3"
        assert before["result"]["status"] == after["result"]["status"] == "ok"
        for key in COLLECTION_PAIR_IDENTITY_KEYS:
            assert before["identity"][key] == after["identity"][key]
        assert before["identity"]["benchmark_sha256"] == (
            PROFILE_BENCHMARK_SHA256 if comparison["allocation_instrumented"] else ORDINARY_BENCHMARK_SHA256)
        for key in COLLECTION_RESULT_KEYS:
            assert before["result"][key] == after["result"][key]
        assert changed_sources(before, after) == {"isovar/__init__.py", "isovar/read_collector.py"}
        modes.add(after["identity"]["mode"])
    assert modes == {"defaults", "primary_only"}
    # These assertions check the frozen observations, not machine-dependent
    # runtime or memory thresholds for code executed during CI.
    for label in ("t3dedup", "t3dedupprimary", "t3tagged"):
        before, after = [data["runs"][phase + "_" + label]["result"] for phase in ("baseline", "final")]
        assert after["locus_collection_peak_rss_bytes"] < 0.35 * before["locus_collection_peak_rss_bytes"]


def test_full_pipeline_ranked_proteins_and_rna_are_unchanged(collection_evidence):
    data = collection_evidence["pipeline"]
    assert len(data["comparisons"]) == 7
    expected_windows = dict(t1primary=2006, t3dedup=2854, t3tagged=4351,
                            ontdync=241, ontdel=384, bulk=112, empty=0)
    for comparison in data["comparisons"]:
        assert comparison["exact_rna_and_ranked_proteins"] is True
        before, after = [data["runs"][comparison[k]] for k in ("before", "after")]
        assert before["identity"]["isovar"] == "1.8.2"
        assert after["identity"]["isovar"] == "1.8.3"
        for run in (before, after):
            assert run["result"]["status"] == run["result"]["validation_status"] == "ok"
            assert all(m["status"] == "ok" for m in run["measurements"].values())
        for key in PIPELINE_RESULT_KEYS:
            assert before["result"][key] == after["result"][key]
        label = comparison["after"].removeprefix("final_")
        assert after["result"]["checked_protein_count"] == expected_windows[label]
    assert sum(expected_windows.values()) == 9948


def test_allocation_profiles_and_failed_runs_are_explicit(collection_evidence):
    runs = collection_evidence["collection"]["runs"]
    before, after = [runs[name]["allocations"] for name in ("baseline_profile_t1", "final_profile_t1")]
    assert before["stage"] == after["stage"] == "before_mate_merging"
    assert before["locus_read_count"] == after["locus_read_count"] == 34071
    assert before["query_bases"] == after["query_bases"] == 30619733
    assert after["traced_current_bytes"] < 0.3 * before["traced_current_bytes"]
    for profile in (before, after):
        assert profile["tracer_overhead_bytes"] > 0
        assert profile["traced_peak_bytes"] >= profile["traced_current_bytes"]
        assert profile["locations"]
    assert runs["interrupted_profile"]["result"]["status"] == "error"
    failed = collection_evidence["pipeline"]["runs"]["interrupted_output"]["result"]
    assert failed["status"] == "error" and failed["counts"] is None
    assert failed["result_record_recovered_from_exception_log"] is True
    assert "No space left on device" in failed["error"]
