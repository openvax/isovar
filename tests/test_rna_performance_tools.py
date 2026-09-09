"""Performance evidence must fail closed on changed inputs or incomplete runs."""

from copy import deepcopy
import json
from types import SimpleNamespace

import pytest

from isovar.allele_read import AlleleRead
from isovar.variant_sequence import VariantSequence
from tests.data.osteosarc.expansion import benchmark, performance_report
from tests.data.osteosarc.expansion.inventory import write_json


def test_rna_fingerprint_preserves_identity_multiplicity_coverage_and_order():
    read = AlleleRead("AAA", "C", "GGG", "same", source_read_count=2)
    sequence = VariantSequence("AA", "C", "GG", [read])
    expected = benchmark.rna_fingerprint([sequence])
    assert benchmark.rna_fingerprint([deepcopy(sequence)]) == expected
    for name, value in [("prefix", "TAAA"), ("suffix", "GGGT"), ("allele", "T"),
                        ("name", "different"), ("source_read_count", 1)]:
        values = {k: getattr(read, k) for k in ("prefix", "allele", "suffix", "name", "source_read_count")}
        values[name] = value
        altered = VariantSequence("AA", "C", "GG", [AlleleRead(**values)])
        assert benchmark.rna_fingerprint([altered]) != expected
    longer = VariantSequence("AAA", "C", "GGG", [read])
    assert benchmark.rna_fingerprint([longer, sequence]) != benchmark.rna_fingerprint([sequence, longer])


@pytest.mark.parametrize("timeout", [0, -1, float("inf"), float("nan"), 1801])
def test_invalid_benchmark_budget_fails_before_creating_output(timeout):
    with pytest.raises(ValueError, match="finite timeout"):
        benchmark.run(SimpleNamespace(timeout=timeout))


def test_measurement_phases_are_separate_and_keep_failures(tmp_path):
    measurements = {}
    with benchmark.measure("production", tmp_path, measurements, 1, False):
        pass
    with pytest.raises(ValueError, match="test error"):
        with benchmark.measure("validation", tmp_path, measurements, 1, False):
            raise ValueError("test error")
    assert json.loads((tmp_path / "production.json").read_text())["status"] == "ok"
    assert json.loads((tmp_path / "validation.json").read_text())["status"] == "interrupted"
    assert set(measurements) == {"production", "validation"}


def benchmark_pair(tmp_path):
    identity = dict(input_bam_sha256="bam", input_index_sha256="index", variant_id="variant",
                    mode="defaults", reference_manifest_sha256="reference")
    result = dict(status="ok", validation_status="ok", rna_sha256="rna", proteins_sha256="proteins",
                  counts={"alt": 2}, public_protein_count=1, checked_protein_count=20)
    directories = [tmp_path / label for label in ("before", "after")]
    for directory in directories:
        directory.mkdir()
        write_json(directory / "identity.json", identity)
        write_json(directory / "result.json", result)
        write_json(directory / "production.json", dict(wall_seconds=1))
    entries = [f"{p.name}={p}" for p in directories]
    return entries, directories


def test_benchmark_comparison_validates_complete_outputs(tmp_path):
    entries, _ = benchmark_pair(tmp_path)
    result = performance_report.collect_benchmarks(entries, ["before=after"])
    assert result["comparisons"] == [dict(before="before", after="after", exact_rna_and_ranked_proteins=True)]


@pytest.mark.parametrize("file,key,value", [
    ("identity", "mode", "primary_only"), ("identity", "input_bam_sha256", "different"),
    ("identity", "reference_manifest_sha256", "different"),
    ("result", "status", "resource_limit"), ("result", "validation_status", "error"),
    ("result", "rna_sha256", "different"), ("result", "proteins_sha256", "different"),
    ("result", "counts", {"alt": 1}), ("result", "checked_protein_count", 1),
])
def test_benchmark_comparison_rejects_changed_or_incomplete_evidence(tmp_path, file, key, value):
    entries, directories = benchmark_pair(tmp_path)
    path = directories[1] / (file + ".json")
    data = json.loads(path.read_text())
    data[key] = value
    path.write_text(json.dumps(data))
    with pytest.raises(ValueError):
        performance_report.collect_benchmarks(entries, ["before=after"])


@pytest.mark.parametrize("n_names", [0, 1, 7, 8, 9, 31, 100])
def test_support_bitsets_round_trip_without_modifying_input(n_names):
    names = [f"read-{i:03d}" for i in range(n_names)]
    proteins = [dict(supporting_read_names=names[::step], checks=[2, 1]) for step in (1, 2, 3)]
    row = dict(defaults=dict(proteins=proteins[:1], uncapped_ranked_proteins=proteins),
               primary_only=dict(proteins=proteins[1:]))
    untouched = deepcopy(row)
    packed = performance_report.pack_support_names(row)
    assert row == untouched
    for mode in ("defaults", "primary_only"):
        for field in row[mode]:
            for original, protein in zip(row[mode][field], packed[mode][field]):
                assert performance_report.supporting_read_names(packed, protein) == original["supporting_read_names"]
                assert protein["checks"] == [1, 2]


@pytest.mark.parametrize("changed,value", [("bits", ""), ("bits", "Ag=="), ("bits", "!invalid!"),
                                          ("names", ["different"]), ("encoding", "unknown")])
def test_support_bitsets_reject_corruption(changed, value):
    row = dict(defaults=dict(proteins=[dict(supporting_read_names=["original"], checks=[])]), primary_only={})
    packed = performance_report.pack_support_names(row)
    protein = packed["defaults"]["proteins"][0]
    if changed == "bits":
        protein["supporting_read_name_bits"] = value
    elif changed == "names":
        packed["supporting_read_name_table"] = value
    else:
        packed["support_encoding"] = value
    with pytest.raises(ValueError):
        performance_report.supporting_read_names(packed, protein)


def rerun_inputs(tmp_path, monkeypatch):
    monkeypatch.setattr(performance_report, "SOURCE_IDS", ["source"])
    audit, baseline = tmp_path / "audit", tmp_path / "baseline"
    for directory in (audit, baseline):
        (directory / "source").mkdir(parents=True)
    write_json(audit / "source/run.json", {})
    protein = dict(supporting_read_names=["original"], checks=[dict(status="ok")])
    mode = dict(status="ok", outcome="expected_top_protein", counts={"alt": 2},
                uncapped_validation_status="ok", proteins=[protein], uncapped_ranked_proteins=[protein])
    for variant_id in [performance_report.MT_VARIANT, *(f"nuclear-{i}" for i in range(43))]:
        current = dict(variant_id=variant_id, source_id="source", software="new",
                       defaults=deepcopy(mode), primary_only=deepcopy(mode), source_bam_sha256="bam")
        previous = deepcopy(current)
        previous["software"] = "old"
        if variant_id == performance_report.MT_VARIANT:
            for m in ("defaults", "primary_only"):
                previous[m].update(status="resource_limit", outcome="audit_timeout")
        write_json(audit / f"source/{variant_id}.json", current)
        write_json(baseline / f"source/{variant_id}.json", previous)
    return audit, baseline


def test_rerun_packaging_keeps_full_support_and_checks_all_controls(tmp_path, monkeypatch):
    result = performance_report.collect_reruns(*rerun_inputs(tmp_path, monkeypatch))
    assert result["unchanged_nuclear_rows"] == 43
    assert result["unchanged_nuclear_modes"] == 86
    assert result["recovered_mitochondrial_modes"] == 2
    assert len(result["checkpoints"]) == 44
    row, = result["rows"]
    assert performance_report.supporting_read_names(row, row["defaults"]["proteins"][0]) == ["original"]


@pytest.mark.parametrize("change", ["counts", "timeout", "validation", "empty", "input", "nuclear"])
def test_rerun_packaging_rejects_missing_or_changed_evidence(tmp_path, monkeypatch, change):
    audit, baseline = rerun_inputs(tmp_path, monkeypatch)
    variant_id = "nuclear-0" if change == "nuclear" else performance_report.MT_VARIANT
    path = audit / f"source/{variant_id}.json"
    row = json.loads(path.read_text())
    if change == "counts":
        row["defaults"]["counts"] = {"alt": 1}
    elif change == "timeout":
        row["defaults"]["status"] = "resource_limit"
    elif change == "validation":
        row["defaults"]["uncapped_validation_status"] = "error"
    elif change == "empty":
        row["defaults"]["proteins"] = []
    elif change == "input":
        row["source_bam_sha256"] = "different"
    else:
        row["defaults"]["proteins"][0]["supporting_read_names"] = ["different"]
    path.write_text(json.dumps(row))
    with pytest.raises(ValueError):
        performance_report.collect_reruns(audit, baseline)
