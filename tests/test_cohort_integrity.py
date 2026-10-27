"""Reports must describe the scientific run whose bytes they actually consume."""

from copy import deepcopy
from dataclasses import asdict
import json
from pathlib import Path
import shutil

import pytest

from isovar.cohort import run_cohort, _additional_artifacts, _check_additional_receipt
from isovar.cli.isovar_cohort import run as cohort_cli
from isovar.sv_rna_audit import file_sha256, read_json, write_json
from .test_fusion import example
from .testing_helpers import data_path


@pytest.fixture(scope="module")
def recorded_small_cohort(tmp_path_factory):
    root = tmp_path_factory.mktemp("recorded-cohort")
    calls = root / "calls.json"
    write_json(calls, [dict(contig="14", start=92811641, ref="GAGA", alt="G")])
    bam = Path(data_path("data/colo829_whole_run/SRR8615617.all_sites.grch38.bam"))
    manifest = dict(schema="isovar.cohort_manifest.v1", samples=[dict(sample_id="COLO829")],
        files={name: dict(path=str(path), sha256=file_sha256(path)) for name, path in
               [("bam", bam), ("calls", calls)]},
        sources=[dict(source_id="original", sample_id="COLO829", library_id="colo-library", bam="bam",
                      scope="Regional RNA", provenance={"accession": "SRR8615617"})],
        small_variants=[dict(source_id="original", variants="calls", ensembl_release=87,
                             protein_settings=dict(protein_sequence_length=33))])
    path = root / "manifest.json"
    write_json(path, manifest)
    run_cohort(path, root / "out", render=False)
    return manifest, root / "out"


@pytest.mark.parametrize("change", ["sample", "library", "input", "reference", "settings", "bam"])
def test_changed_manifest_cannot_relabel_recorded_small_variants(recorded_small_cohort, tmp_path, change):
    original, output = recorded_small_cohort
    manifest = deepcopy(original)
    shutil.copytree(output, tmp_path / "out")
    if change == "sample":
        manifest["samples"][0]["sample_id"] = "HCC1395"
        manifest["sources"][0]["sample_id"] = "HCC1395"
    elif change == "library":
        manifest["sources"][0]["library_id"] = "another-library"
    elif change in ("input", "bam"):
        manifest["files"]["calls" if change == "input" else "bam"]["sha256"] = "0" * 64
    elif change == "reference":
        manifest["small_variants"][0]["ensembl_release"] = 75
    else:
        manifest["small_variants"][0]["protein_settings"]["protein_sequence_length"] = 99
    path = tmp_path / "manifest.json"
    write_json(path, manifest)
    saved = {p.name: p.read_bytes() for p in (tmp_path / "out").glob("*.json")}
    with pytest.raises(ValueError, match="analysis identity differs"):
        run_cohort(path, tmp_path / "out", render_only=True, render=False)
    assert saved == {p.name: p.read_bytes() for p in (tmp_path / "out").glob("*.json")}


def test_render_allows_relocation_and_presentation_without_original_bam(recorded_small_cohort, tmp_path):
    original, output = recorded_small_cohort
    manifest = deepcopy(original)
    shutil.copytree(output, tmp_path / "out")
    for value in manifest["files"].values():
        value["path"] = "offline-input"
    manifest["small_variants"][0]["title"] = "A new display title"
    manifest["samples"][0]["title"] = "New sample display title"
    path = tmp_path / "manifest.json"
    write_json(path, manifest)
    assert run_cohort(path, tmp_path / "out", render_only=True, render=False) == []
    summary = read_json(tmp_path / "out/cohort.json")
    assert summary["small_variants"][0]["identity"]["source"]["sample_id"] == "COLO829"


def test_cli_reports_identity_mismatch_and_corrupt_small_result(recorded_small_cohort, tmp_path, capsys):
    original, output = recorded_small_cohort
    manifest = deepcopy(original)
    shutil.copytree(output, tmp_path / "out")
    path = tmp_path / "manifest.json"
    manifest["sources"][0]["library_id"] = "changed"
    write_json(path, manifest)
    with pytest.raises(SystemExit) as error:
        cohort_cli(["--manifest", str(path), "--output-dir", str(tmp_path / "out"), "--render-only"])
    assert error.value.code == 2
    assert "analysis identity differs" in capsys.readouterr().err
    write_json(path, original)
    (tmp_path / "out/sources/original/default/variants.jsonl.gz").write_bytes(b"corrupt")
    with pytest.raises(ValueError, match="result digest mismatch"):
        run_cohort(path, tmp_path / "out", render_only=True, render=False)


def fusion_manifest(root, **example_settings):
    fusion, refs, reads = example(**example_settings)
    supplied = dict(fusion=asdict(fusion), references=[asdict(r) for r in refs], reads=[asdict(r) for r in reads])
    path = root / "fusion.json"
    write_json(path, supplied)
    return dict(schema="isovar.cohort_manifest.v1", samples=[dict(sample_id="sample")],
        files={"fusion": dict(path=str(path), sha256=file_sha256(path))},
        sources=[dict(source_id="original", sample_id="sample", library_id="library",
                      scope="Supplied fusion", provenance={"source": "synthetic regression"})],
        fusions=[dict(id="fusion", sample_id="sample", source_id="original", input="fusion")])


@pytest.mark.parametrize("change", ["sample", "library", "input", "settings", "result", "missing-receipt"])
def test_additional_receipts_reject_changed_identity_or_bytes(tmp_path, change):
    manifest = fusion_manifest(tmp_path)
    path, output = tmp_path / "manifest.json", tmp_path / "out"
    write_json(path, manifest)
    run_cohort(path, output, render=False)
    if change == "sample":
        manifest["samples"][0]["sample_id"] = manifest["sources"][0]["sample_id"] = "different"
        manifest["fusions"][0]["sample_id"] = "different"
    elif change == "library":
        manifest["sources"][0]["library_id"] = "different"
    elif change == "input":
        manifest["files"]["fusion"]["sha256"] = "0" * 64
    elif change == "settings":
        manifest["fusions"][0]["parameters"] = dict(min_fragments=100)
    elif change == "result":
        write_json(output / "fusions/fusion/result.json.gz", {"modified": True})
    else:
        (output / "fusions/fusion/receipt.json").unlink()
    write_json(path, manifest)
    with pytest.raises(ValueError, match="identity differs|digest mismatch|Missing additional analysis receipt"):
        run_cohort(path, output, render_only=True, render=False)


def test_additional_resume_presentation_and_input_sample_validation(tmp_path):
    manifest = fusion_manifest(tmp_path)
    path, output = tmp_path / "manifest.json", tmp_path / "out"
    write_json(path, manifest)
    run_cohort(path, output, render=False)
    result = output / "fusions/fusion/result.json.gz"
    before = result.stat().st_mtime_ns
    run_cohort(path, output, resume=True, render=False)
    assert result.stat().st_mtime_ns == before
    manifest["fusions"][0]["title"] = "Revised display title"
    write_json(path, manifest)
    run_cohort(path, output, render_only=True, render=False)
    manifest["samples"][0]["sample_id"] = manifest["sources"][0]["sample_id"] = "different"
    manifest["fusions"][0]["sample_id"] = "different"
    write_json(path, manifest)
    with pytest.raises(ValueError, match="fusion and cohort source"):
        run_cohort(path, tmp_path / "new-output", render=False)


def test_additional_tasks_require_sample_in_api_and_cli(tmp_path):
    manifest = fusion_manifest(tmp_path)
    del manifest["fusions"][0]["sample_id"]
    path = tmp_path / "manifest.json"
    write_json(path, manifest)
    with pytest.raises(ValueError, match="require sample_id"):
        run_cohort(path, tmp_path / "out", render=False)


def test_sv_receipt_includes_referenced_path_artifacts(tmp_path):
    target = tmp_path / "report.json"
    write_json(target, dict(runs=[dict(artifact="path.json.gz")]))
    write_json(tmp_path / "path.json.gz", dict(sequence="ACGT"))
    identity = dict(task={}, source={}, inputs={})
    receipt = dict(schema="isovar.cohort_additional.v1", identity=identity,
                   artifacts=_additional_artifacts("sv_audits", tmp_path, target))
    path = tmp_path / "receipt.json"
    write_json(path, receipt)
    assert _check_additional_receipt(path, identity, tmp_path, target)
    write_json(tmp_path / "path.json.gz", dict(sequence="AAAA"))
    with pytest.raises(ValueError, match="artifact digest mismatch"):
        _check_additional_receipt(path, identity, tmp_path, target)


def test_prior_assessment_snapshot_is_verified_without_original_input(tmp_path):
    manifest = fusion_manifest(tmp_path)
    manifest["prior_assessments"] = manifest.pop("fusions")
    path, output = tmp_path / "manifest.json", tmp_path / "out"
    write_json(path, manifest)
    run_cohort(path, output, render=False)
    Path(manifest["files"]["fusion"]["path"]).unlink()
    assert run_cohort(path, output, render_only=True, render=False) == []
    summary = json.loads((output / "cohort.json").read_text())
    assert summary["additional"][0]["status"] == "recorded"
