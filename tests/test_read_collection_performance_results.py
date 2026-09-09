"""Audit the recorded full-depth #229 evidence without loading huge BAMs."""

import json
from pathlib import Path

import pytest

from tests.data.osteosarc.expansion.inventory import digest


RESULTS = Path(__file__).parent / "data/osteosarc/expansion/performance/collection-results"


@pytest.fixture(scope="module")
def collection_evidence():
    manifest = json.loads((RESULTS / "manifest.json").read_text())
    for name, checksum in manifest["files"].items():
        assert digest(RESULTS / name) == checksum
    return {name: json.loads((RESULTS / (name + ".json")).read_text()) for name in ("collection", "pipeline")}


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
        for key in ("input_bam_sha256", "input_index_sha256", "inventory_sha256", "variant_id", "mode",
                    "merge_overlapping_fragments", "allocation_instrumented"):
            assert before["identity"][key] == after["identity"][key]
        for key in ("locus_reads_sha256", "allele_reads_sha256", "locus_read_count", "allele_read_count", "counts"):
            assert before["result"][key] == after["result"][key]
        changes = {p for p, sha in before["identity"]["source_files"].items()
                   if after["identity"]["source_files"][p] != sha}
        assert changes == {"isovar/__init__.py", "isovar/read_collector.py"}
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
        for key in ("rna_sha256", "proteins_sha256", "counts", "public_protein_count", "checked_protein_count"):
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
