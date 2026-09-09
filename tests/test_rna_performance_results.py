"""Offline checks of the full-depth #227 evidence, not reruns of giant BAMs."""

import gzip
import json
from pathlib import Path

import pytest

from tests.data.osteosarc.expansion.inventory import digest
from tests.data.osteosarc.expansion.performance_report import SOURCE_IDS, supporting_read_names


RESULTS = Path(__file__).parent / "data/osteosarc/expansion/performance/results"


@pytest.fixture(scope="module")
def performance_evidence():
    manifest = json.loads((RESULTS / "manifest.json").read_text())
    for name, checksum in manifest["files"].items():
        assert digest(RESULTS / name) == checksum
    with gzip.open(RESULTS / "reruns.json.gz", "rt") as handle:
        return json.load(handle)


def test_full_depth_reruns_recover_every_interrupted_mode(performance_evidence):
    data = performance_evidence
    assert data["recovered_mitochondrial_modes"] == 10
    assert data["unchanged_nuclear_rows"] == 215
    assert data["unchanged_nuclear_modes"] == 430
    assert len(data["checkpoints"]) == 220
    assert {row["source_id"] for row in data["rows"]} == set(SOURCE_IDS)
    windows, checks = 0, 0
    for row in data["rows"]:
        assert row["variant_id"] == "MT_ND5-chrM-12994"
        assert row["mitochondrial_interpretation"]["numt_origin_status"] == "not_ruled_out"
        for mode in ("defaults", "primary_only"):
            result = row[mode]
            assert result["status"] == result["uncapped_validation_status"] == "ok"
            assert result["outcome"] == "expected_top_protein"
            ranked = result["uncapped_ranked_proteins"]
            assert result["returned_protein_limit"] == 1
            assert result["proteins"] == ranked[:1]
            assert 26 <= len(ranked[0]["amino_acids"]) <= 28
            assert max(len(p["amino_acids"]) for p in ranked) == 49
            windows += len(ranked)
            for protein in ranked:
                assert protein["validation_status"] == "ok"
                assert protein["checks"]
                assert all(c["status"] == "ok" and c["genetic_code"] == 2 for c in protein["checks"])
                checks += len(protein["checks"])
            # Decode edge and middle ranks, preserving exact names (not just
            # trusting a count). The full payload is checksum-pinned above.
            for protein in (ranked[0], ranked[len(ranked) // 2], ranked[-1]):
                assert len(supporting_read_names(row, protein)) == protein["supporting_template_names"]
    assert windows == 28890
    assert checks == 64195


def test_audit_source_and_settings_provenance_is_consistent(performance_evidence):
    runs = list(performance_evidence["runs"].values())
    for run in runs:
        assert run["software"]["isovar"] == "1.8.2"
        assert run["time_limit_seconds_per_mode"] == 600
        for field in ("sources", "settings", "inventory_sha256", "reference_manifest_sha256"):
            assert run[field] == runs[0][field]
    benchmarks = json.loads((RESULTS / "benchmarks.json").read_text())["runs"]
    for path, checksum in benchmarks["final_t3"]["identity"]["source_files"].items():
        if path in runs[0]["sources"]:
            assert checksum == runs[0]["sources"][path]


def test_completed_full_depth_comparisons_preserve_fingerprints():
    data = json.loads((RESULTS / "benchmarks.json").read_text())
    assert len(data["comparisons"]) == 5
    for comparison in data["comparisons"]:
        assert comparison["exact_rna_and_ranked_proteins"] is True
        before, after = [data["runs"][comparison[key]] for key in ("before", "after")]
        assert before["identity"]["isovar"] == "1.8.1"
        assert after["identity"]["isovar"] == "1.8.2"
        for key in ("rna_sha256", "proteins_sha256", "counts", "checked_protein_count"):
            assert before["result"][key] == after["result"][key]
    # Retain failed intermediate measurements rather than silently dropping
    # them from the evidence or calling their counts zero.
    failed = data["runs"]["interval_coverage_timeout"]["result"]
    assert failed["status"] == "resource_limit"
    assert failed["counts"]["reads"]["alt"] == 29130
