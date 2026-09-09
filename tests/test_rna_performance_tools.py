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
