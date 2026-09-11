"""Audit the recorded full-depth #234 coordinate-sharing evidence without huge BAMs."""

import json
from pathlib import Path

import pytest

from tests.data.osteosarc.expansion.inventory import digest


EXPANSION = Path(__file__).parent / "data/osteosarc/expansion"
RESULTS = EXPANSION / "performance/coordinate-sharing-results"
COLLECTION_LABELS = {"t1", "t1primary", "t3dedup", "t3dedupprimary", "t3tagged", "ontdync", "ontdel", "bulk",
                     "t3exoc4", "profile_t1"}


@pytest.fixture(scope="module")
def evidence():
    manifest = json.loads((RESULTS / "manifest.json").read_text())
    assert set(manifest) == {"files", "generator_sha256"}
    for name, checksum in manifest["files"].items():
        assert digest(RESULTS / name) == checksum
    # The shipped validator must be the one that reproduces this package.
    assert manifest["generator_sha256"] == digest(EXPANSION / "collection_report.py")
    return {name: json.loads((RESULTS / (name + ".json")).read_text()) for name in ("collection", "pipeline")}


def changed_isovar_sources(before, after):
    return {path for path, checksum in before["identity"]["source_files"].items()
            if path.startswith("isovar/") and after["identity"]["source_files"].get(path) != checksum}


def test_per_call_sharing_preserves_every_full_depth_collection(evidence):
    data = evidence["collection"]
    assert {c["after"].removeprefix("final_") for c in data["comparisons"]} == COLLECTION_LABELS
    for comparison in data["comparisons"]:
        assert comparison["exact_locus_and_allele_reads"] is True
        before, after = (data["runs"][comparison[key]] for key in ("before", "after"))
        assert (before["identity"]["isovar"], after["identity"]["isovar"]) == ("1.8.3", "1.8.4")
        assert before["identity"]["benchmark_sha256"] == after["identity"]["benchmark_sha256"]
        assert changed_isovar_sources(before, after) == {"isovar/__init__.py", "isovar/read_collector.py"}
        for key in ("locus_reads_sha256", "allele_reads_sha256", "locus_read_count", "allele_read_count", "counts"):
            assert before["result"][key] == after["result"][key]


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
        for run in (before, after):
            assert run["result"]["status"] == run["result"]["validation_status"] == "ok"
            assert all(m["status"] == "ok" for m in run["measurements"].values())
        for key in ("rna_sha256", "proteins_sha256", "counts", "public_protein_count", "checked_protein_count"):
            assert before["result"][key] == after["result"][key]
        assert after["result"]["checked_protein_count"] == 4351


def test_frozen_collection_observations_keep_memory_and_lower_total_cpu(evidence):
    # These check the recorded observations, not machine-dependent thresholds
    # for code executed during CI.
    runs = evidence["collection"]["runs"]
    totals = dict(baseline=0.0, final=0.0)
    for label in COLLECTION_LABELS - {"profile_t1"}:
        before, after = (runs[f"{phase}_{label}"]["result"] for phase in ("baseline", "final"))
        assert abs(after["locus_collection_peak_rss_bytes"] / before["locus_collection_peak_rss_bytes"] - 1) < 0.05
        totals["baseline"] += before["locus_collection_cpu_seconds"]
        totals["final"] += after["locus_collection_cpu_seconds"]
    assert totals["final"] < 0.95 * totals["baseline"]


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
