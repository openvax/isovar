"""The reporting harness must distinguish errors, zeros and analysis stages."""

from collections import Counter
import gzip
import json
from pathlib import Path
from types import SimpleNamespace

import pysam
import pytest

from .data.osteosarc import audit


DATA = Path(__file__).parent / "data" / "osteosarc"
MANIFEST = json.loads((DATA / "manifest.json").read_text())


def test_audit_collection_exception_has_null_counts_not_biological_zero(monkeypatch):
    def fail(*args, **kwargs):
        raise RuntimeError("deliberate collection failure")

    monkeypatch.setattr(audit.ReadCollector, "read_evidence_for_variant", fail)
    result = audit.isovar_result(None, None)
    assert result == {"status": "error", "stage": "read_evidence", "counts": None,
                      "error": {"type": "RuntimeError", "message": "deliberate collection failure"}}
    assert "sequence_creation_default" not in result


def test_audit_sequence_exception_preserves_successful_allele_counts(monkeypatch):
    evidence = SimpleNamespace(ref_reads=[1, 2], alt_reads=[3], other_reads=[], alt_read_names={"alt"})
    monkeypatch.setattr(audit.ReadCollector, "read_evidence_for_variant", lambda *a: evidence)

    def fail(*args, **kwargs):
        raise ValueError("deliberate sequence failure")

    monkeypatch.setattr(audit.VariantSequenceCreator, "reads_to_variant_sequences", fail)
    result = audit.isovar_result(None, None)
    assert result["status"] == "ok"
    assert result["counts"] == {"ref": 2, "alt": 1, "other": 0}
    for key in ("sequence_creation_default", "sequence_creation_coverage1"):
        assert result[key]["status"] == "error"
        assert result[key]["stage"] == "sequence_creation"
        assert result[key]["error"]["type"] == "ValueError"


def test_audit_genuine_empty_collection_is_success_not_exception(monkeypatch):
    evidence = SimpleNamespace(ref_reads=[], alt_reads=[], other_reads=[], alt_read_names=set())
    monkeypatch.setattr(audit.ReadCollector, "read_evidence_for_variant", lambda *a: evidence)
    result = audit.isovar_result(None, None)
    assert result["status"] == "ok"
    assert result["counts"] == {"ref": 0, "alt": 0, "other": 0}
    assert result["sequence_creation_default"]["candidate_count"] == 0


def test_audit_rejects_source_drift_before_producing_counts(tmp_path):
    (tmp_path / "bulk_star_t0.sam").write_text("not the pinned source\n")
    with pytest.raises(ValueError, match="checksum mismatch"):
        audit.build_report(tmp_path)


@pytest.fixture(scope="module")
def fixture_results(tmp_path_factory):
    directory = tmp_path_factory.mktemp("osteosarc_audit")
    result = {}
    for name, data in MANIFEST["datasets"].items():
        sam = directory / (name + ".sam")
        sam.write_bytes(gzip.decompress((DATA / data["file"]).read_bytes()))
        full = directory / (name + ".bam")
        primary = directory / (name + ".primary.bam")
        pysam.sort("--no-PG", "-o", str(full), str(sam))
        pysam.index(str(full))
        pysam.view("--no-PG", "-b", "-F", str(audit.PRIMARY_EXCLUDE_FLAGS), "-o",
                   str(primary), str(full), catch_stdout=False)
        pysam.index(str(primary))
        result[name] = audit.audit_alignment(name, full, primary, data["variants"])
    return result


@pytest.mark.parametrize("name", MANIFEST["datasets"])
def test_audit_completes_other_loci_when_real_spliced_deletions_crash(name, fixture_results):
    rows = fixture_results[name]
    assert len(rows) == 6
    for row in rows:
        crashing = name == "ont_t1" and row["variant"]["gene"] in ("PIP5K1A", "H1-2", "GTF3C5")
        for mode in ("primary_only", "defaults"):
            result = row[mode]
            assert result["status"] == ("error" if crashing else "ok")
            if crashing:
                # This pins the report of #217, not a claim of correct Isovar
                # handling. Update this expectation when #217 is fixed.
                assert result["counts"] is None
                assert result["stage"] == "read_evidence"
                assert result["error"] == {"type": "AttributeError", "message": "'NoneType' object has no attribute 'upper'"}


@pytest.mark.parametrize("name", MANIFEST["datasets"])
def test_audit_primary_snp_counts_match_frozen_independent_observations(name, fixture_results):
    for row in fixture_results[name]:
        record = row["variant"]
        if len(record["ref"]) != 1:
            continue
        observations = [o for o in record["observations"]
                        if not o["flag"] & audit.PRIMARY_EXCLUDE_FLAGS and o["mapq"] >= 1]
        counts = Counter(o["allele"] for o in observations)
        expected = {"ref": counts[record["ref"]], "alt": counts[record["alt"]],
                    "other": len(observations) - counts[record["ref"]] - counts[record["alt"]]}
        assert row["primary_only"]["counts"] == row["independent"]["counts"] == expected


def test_full_region_snapshot_is_complete_and_separate_from_selected_fixtures():
    report = json.loads((DATA / "sample_audit.json").read_text())
    expected = {(n, v["variant_id"]) for n, d in MANIFEST["datasets"].items() for v in d["variants"]}
    assert len(report["rows"]) == len(expected) == 12
    assert {(r["sample"], r["variant"]["variant_id"]) for r in report["rows"]} == expected
    for name, source in report["sources"].items():
        assert source["sha256"] == MANIFEST["datasets"][name]["source_sam_sha256"]
        assert source["region_records"] > MANIFEST["datasets"][name]["fixture_records"]
    for row in report["rows"]:
        for mode in ("primary_only", "defaults"):
            result = row[mode]
            if result["status"] == "error":
                assert result["counts"] is None
            else:
                assert all(isinstance(n, int) and n >= 0 for n in result["counts"].values())
                assert result["alt_read_names"] <= result["counts"]["alt"]
                assert result["sequence_creation_default"]["represented_read_names"] <= result["alt_read_names"]
