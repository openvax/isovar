"""Portable evidence must replay without original paths or silent drift."""

from copy import deepcopy
import json
from pathlib import Path
import shutil
import tarfile

import pytest

from isovar.audit_bundle import AuditBundle, AuditBundleError, digest, identity
from isovar.audit_replay import replay_sv_rna, read_json
from isovar.sid_audit_bundle import bundle_sid_audit
from isovar.cli.isovar_audit_bundle import run


def source_file(tmp_path):
    path = tmp_path / "original.txt"
    path.write_bytes(b"original sequence evidence\n")
    return path


def test_bundle_archive_roundtrip_survives_removal_of_originals(tmp_path):
    original = source_file(tmp_path)
    bundle = AuditBundle.create(tmp_path / "bundle", {"inputs/reads.txt": original},
                                metadata={"sample": "T2"}, expected_sha256={"inputs/reads.txt": digest(original)})
    bundle.pack(tmp_path / "backup.tar.gz")
    identity_before = bundle.bundle_id
    original.unlink()
    shutil.rmtree(bundle.directory)
    restored = AuditBundle.unpack(tmp_path / "backup.tar.gz", tmp_path / "relocated")
    assert restored.bundle_id == identity_before
    assert restored.verify() == dict(bundle_id=identity_before, files=1, bytes=27)
    assert restored.path("inputs/reads.txt").read_bytes() == b"original sequence evidence\n"
    manifest = restored.manifest
    manifest["files"].clear()
    assert len(restored.manifest["files"]) == 1


@pytest.mark.parametrize("name", ["../outside", "/absolute", "a/../../outside", "a//b", "./a", "a\\b", ""])
def test_create_rejects_nonportable_file_names_before_writing(tmp_path, name):
    with pytest.raises(AuditBundleError):
        AuditBundle.create(tmp_path / "bundle", {name: source_file(tmp_path)})
    assert not (tmp_path / "bundle").exists()


def test_bundle_does_not_replace_an_existing_destination(tmp_path):
    destination = tmp_path / "bundle"
    destination.mkdir()
    sentinel = destination / "user-data"
    sentinel.write_text("preserve")
    with pytest.raises(FileExistsError):
        AuditBundle.create(destination, {"input": source_file(tmp_path)})
    assert sentinel.read_text() == "preserve"


def test_interrupted_copy_publishes_no_complete_or_partial_bundle(tmp_path, monkeypatch):
    def interrupted(*args):
        raise KeyboardInterrupt()
    monkeypatch.setattr(shutil, "copyfile", interrupted)
    with pytest.raises(KeyboardInterrupt):
        AuditBundle.create(tmp_path / "bundle", {"input": source_file(tmp_path)})
    assert not (tmp_path / "bundle").exists()


def test_expected_input_pin_is_enforced_during_publication(tmp_path):
    with pytest.raises(AuditBundleError, match="Input checksum mismatch"):
        AuditBundle.create(tmp_path / "bundle", {"input": source_file(tmp_path)},
                            expected_sha256={"input": "0" * 64})
    assert not (tmp_path / "bundle").exists()


def test_opened_bundle_detects_an_on_disk_manifest_change(tmp_path):
    bundle = AuditBundle.create(tmp_path / "bundle", {"input": source_file(tmp_path)})
    manifest = bundle.manifest
    manifest["metadata"]["sample"] = "different"
    (bundle.directory / "manifest.json").write_text(json.dumps(manifest))
    with pytest.raises(AuditBundleError, match="manifest changed"):
        bundle.pack(tmp_path / "backup.tar.gz")
    assert not (tmp_path / "backup.tar.gz").exists()


@pytest.mark.parametrize("change", ["same_size", "missing", "escape"])
def test_open_rejects_changed_missing_or_escaping_inputs(tmp_path, change):
    original = source_file(tmp_path)
    bundle = AuditBundle.create(tmp_path / "bundle", {"input": original})
    copied = bundle.path("input")
    copied.unlink()
    if change == "same_size":
        copied.write_bytes(b"X" * original.stat().st_size)
    elif change == "escape":
        copied.symlink_to(original)
    with pytest.raises(AuditBundleError):
        AuditBundle.open(bundle.directory)


@pytest.mark.parametrize("change", ["metadata", "schema"])
def test_open_rejects_changed_manifest_or_unsupported_schema(tmp_path, change):
    bundle = AuditBundle.create(tmp_path / "bundle", {"input": source_file(tmp_path)})
    manifest = bundle.manifest
    if change == "metadata":
        manifest["metadata"]["sample"] = "changed"
    else:
        manifest["schema_version"] = 99
    (bundle.directory / "manifest.json").write_text(json.dumps(manifest))
    with pytest.raises(AuditBundleError):
        AuditBundle.open(bundle.directory)


@pytest.mark.parametrize("change", ["missing", "extra", "duplicate", "link", "changed_bytes"])
def test_unpack_rejects_incomplete_or_invalid_archives(tmp_path, change):
    bundle = AuditBundle.create(tmp_path / "bundle", {"input": source_file(tmp_path)})
    archive = tmp_path / "bad.tar.gz"
    with tarfile.open(archive, "w:gz") as out:
        out.add(bundle.directory / "manifest.json", arcname="manifest.json")
        if change != "missing":
            if change == "changed_bytes":
                bundle.path("input").write_bytes(b"X" * bundle.path("input").stat().st_size)
            out.add(bundle.path("input"), arcname="files/input")
        if change in ("extra", "duplicate"):
            out.add(bundle.path("input"), arcname="../outside" if change == "extra" else "files/input")
        if change == "link":
            link = tarfile.TarInfo("files/link")
            link.type = tarfile.SYMTYPE
            link.linkname = "../../outside"
            out.addfile(link)
    with pytest.raises(AuditBundleError):
        AuditBundle.unpack(archive, tmp_path / "restore")
    assert not (tmp_path / "restore").exists()
    assert not (tmp_path / "outside").exists()


@pytest.fixture
def original_sid_audit(tmp_path):
    # Real, checksum-pinned original SAM records, independent of bundle code.
    from tests.data.fusions.audit_three_fusions import DATA, run_case, pinned_json
    from examples.sid_sv_audit.inventory import write_json

    directory = tmp_path / "old-audit"
    directory.mkdir()
    gid, sid, orientation = "PARD3B--CDKN2B-AS1-CDKN2B", "ONT-T1", "reverse"
    result, _ = run_case(DATA, gid, sid, orientation, directory)
    source = dict(id=sid, url=result["source"], claims=[])
    inventory = dict(sources={sid: source}, geometries={gid: dict(nominations=[gid])})
    write_json(directory / "inventory.json.gz", inventory)
    inventory_sha = digest(directory / "inventory.json.gz")
    write_json(directory / "inventory-pin.json", dict(sha256=inventory_sha))
    pinned = json.loads((DATA / "manifest.json").read_text())
    models = pinned_json(DATA, pinned["references"])[gid]
    references = dict(models={m["transcript_id"]: m for m in models}, inventory_sha256=inventory_sha)
    write_json(directory / "references.json.gz", references)
    reference_sha = digest(directory / "references.json.gz")
    write_json(directory / "references-pin.json", dict(sha256=reference_sha))
    batch = dict(id="batch", geometries=[gid])
    write_json(directory / "batches.json", [batch])
    provenance = directory / "query-names.txt"
    provenance.write_text("retained upstream query-name asset\n")
    acquired = dict(source_id=sid, batch_id="batch", status="incomplete", receipt=dict(limits=["max_rounds"]),
        request=dict(batch=batch, source_identity=identity(source)),
        bam=dict(path=str(directory / "reads.bam"), sha256=digest(directory / "reads.bam")),
        index=dict(path=str(directory / "reads.bam.bai"), sha256=digest(directory / "reads.bam.bai")),
        upstream_receipt_files=[dict(path=str(provenance), sha256=digest(provenance))])
    write_json(directory / "sources" / sid / "batch.json", acquired)
    flags = ["input_acquisition_incomplete", "input_acquisition_limit:max_rounds"]
    result["limitations"] = sorted(set(result["limitations"]) | set(flags))
    result["upstream_acquisition"] = dict(status="incomplete", receipt_sha256=identity(acquired),
                                           limits=["max_rounds"], failed_queries=[])
    details = "reconstructions/" + sid + "/" + gid + ".reverse.json.gz"
    write_json(directory / details, result)
    request = dict(source_id=sid, geometry_id=gid, orientation=orientation,
                    acquisition_sha256=identity(acquired), inventory_sha256=inventory_sha,
                    references_sha256=reference_sha, isovar_version="original-fixture")
    checkpoint = dict(request=request, status=result["status"], acquisition_status="incomplete",
        input_limitations=flags, reconstruction=dict(path=details, sha256=digest(directory / details)))
    write_json(directory / "results" / sid / (gid + ".reverse.json.gz"), checkpoint)
    forward = {key: value for key, value in checkpoint.items() if key != "reconstruction"}
    forward.update(request=dict(request, orientation="forward"), status="not_assessable", reason="unavailable_fixture_view")
    forward_path = directory / "results" / sid / (gid + ".forward.json.gz")
    write_json(forward_path, forward)
    # A report-local file pin has a different base than audit-relative pins.
    report_table = directory / "report" / "input-files.tsv"
    report_table.parent.mkdir()
    report_table.write_text("path\tsha256\n")
    write_json(directory / "report" / "priority-subsets.json.gz", dict(
        input_files_manifest=dict(path="input-files.tsv", sha256=digest(report_table))))
    return directory, "/".join((sid, gid, orientation)), read_json(directory / details)


def test_original_read_replay_after_relocation_is_exact_and_resumable(tmp_path, original_sid_audit, monkeypatch):
    from isovar import audit_replay
    from isovar.sv_rna_orf_export import export_sv_rna_orfs

    original, run_id, expected = original_sid_audit
    bundle = bundle_sid_audit(original, tmp_path / "bundle")
    bundle.pack(tmp_path / "backup.tar.gz")
    shutil.rmtree(original)
    shutil.rmtree(bundle.directory)
    relocated = AuditBundle.unpack(tmp_path / "backup.tar.gz", tmp_path / "elsewhere")
    checkpoint = replay_sv_rna(relocated, run_id, tmp_path / "replay")
    assert checkpoint["result"] == expected
    assert checkpoint["comparison"] == "matched"
    candidate, = [c for c in export_sv_rna_orfs(checkpoint["result"])["candidates"]
                   if c["amino_acids"] == "MNISNIHISTQKKKKKKSRF"]
    assert candidate["rna_support"]["fragments"] == 4
    assert "input_acquisition_incomplete" in checkpoint["result"]["limitations"]
    assert len(relocated.manifest["metadata"]["outcomes"]) == 2
    assert len(relocated.manifest["metadata"]["sv_rna_runs"]) == 1
    assert relocated.path("audit/report/input-files.tsv").read_text() == "path\tsha256\n"

    def must_not_recompute(*args, **kwargs):
        raise AssertionError("Matching checkpoint should have been reused")
    monkeypatch.setattr(audit_replay, "reconstruct_sv_rna", must_not_recompute)
    assert replay_sv_rna(relocated, run_id, tmp_path / "replay") == checkpoint


def test_sid_import_rejects_changed_original_provenance(tmp_path, original_sid_audit):
    original, _, _ = original_sid_audit
    (original / "query-names.txt").write_text("changed original names")
    with pytest.raises(AuditBundleError, match="checksum mismatch"):
        bundle_sid_audit(original, tmp_path / "bundle")
    assert not (tmp_path / "bundle").exists()


@pytest.mark.parametrize("conflicting", [False, True])
def test_sid_import_deduplicates_hashing_without_accepting_conflicting_pins(
        tmp_path, original_sid_audit, monkeypatch, conflicting):
    from isovar import sid_audit_bundle
    from examples.sid_sv_audit.inventory import write_json

    original, _, _ = original_sid_audit
    shared = original / "query-names.txt"
    pins = [dict(path=str(shared), sha256=digest(shared)) for _ in range(10)]
    if conflicting:
        pins[-1]["sha256"] = "0" * 64
    write_json(original / "repeated-annotation-pins.json", dict(files=pins))
    calls = []
    native_digest = sid_audit_bundle.digest
    def counted(path):
        if Path(path) == shared:
            calls.append(path)
        return native_digest(path)
    monkeypatch.setattr(sid_audit_bundle, "digest", counted)
    if conflicting:
        with pytest.raises(AuditBundleError, match="checksum mismatch"):
            bundle_sid_audit(original, tmp_path / "bundle")
        assert not (tmp_path / "bundle").exists()
    else:
        bundle = bundle_sid_audit(original, tmp_path / "bundle")
        assert bundle.path("audit/query-names.txt").read_bytes() == shared.read_bytes()
    assert len(calls) == 1


def test_replay_refuses_corrupt_checkpoint_and_software_drift(tmp_path, original_sid_audit, monkeypatch):
    from isovar import audit_replay

    original, run_id, _ = original_sid_audit
    bundle = bundle_sid_audit(original, tmp_path / "bundle")
    checkpoint = replay_sv_rna(bundle, run_id, tmp_path / "replay")
    path, = (tmp_path / "replay").glob("*.json")
    changed = deepcopy(checkpoint)
    changed["result"]["status"] = "changed"
    path.write_text(json.dumps(changed))
    with pytest.raises(AuditBundleError, match="checksum mismatch"):
        replay_sv_rna(bundle, run_id, tmp_path / "replay")
    monkeypatch.setattr(audit_replay, "software_identity", lambda: {"different": "engine"})
    with pytest.raises(AuditBundleError, match="software differs"):
        replay_sv_rna(bundle, run_id, tmp_path / "new-run")
    with pytest.raises(AuditBundleError, match="request drift"):
        replay_sv_rna(bundle, run_id, tmp_path / "replay", allow_software_change=True)
    comparison = replay_sv_rna(bundle, run_id, tmp_path / "new-run", allow_software_change=True)
    assert comparison["comparison"] == "matched"


def test_replay_mismatch_preserves_full_diagnostic_result(tmp_path, original_sid_audit, monkeypatch):
    from isovar import audit_replay

    original, run_id, result = original_sid_audit
    bundle = bundle_sid_audit(original, tmp_path / "bundle")
    changed = deepcopy(result)
    changed["status"] = "different"
    monkeypatch.setattr(audit_replay, "reconstruct_sv_rna", lambda *args, **kwargs: changed)
    with pytest.raises(AuditBundleError, match="differs from the frozen result"):
        replay_sv_rna(bundle, run_id, tmp_path / "replay")
    path, = (tmp_path / "replay").glob("*.json")
    diagnostic = read_json(path)
    assert diagnostic["comparison"] == "mismatch"
    assert diagnostic["result"]["status"] == "different"


def test_replay_rejects_unknown_runs_and_output_inside_bundle(tmp_path, original_sid_audit):
    original, run_id, _ = original_sid_audit
    bundle = bundle_sid_audit(original, tmp_path / "bundle")
    with pytest.raises(AuditBundleError, match="Unknown SV RNA replay run"):
        replay_sv_rna(bundle, "unknown", tmp_path / "replay")
    with pytest.raises(AuditBundleError, match="outside the immutable bundle"):
        replay_sv_rna(bundle, run_id, bundle.directory / "replay")
    assert not (bundle.directory / "replay").exists()


def test_cli_create_resolves_spec_relative_paths_and_reports_errors(tmp_path, capsys, monkeypatch):
    source_file(tmp_path)
    spec = tmp_path / "inputs.json"
    spec.write_text(json.dumps(dict(files={"input": "original.txt"})))
    elsewhere = tmp_path / "elsewhere"
    elsewhere.mkdir()
    monkeypatch.chdir(elsewhere)
    run(["create", "--spec", str(spec), "--output", str(tmp_path / "bundle")])
    assert json.loads(capsys.readouterr().out)["files"] == 1
    run(["verify", str(tmp_path / "bundle")])
    assert json.loads(capsys.readouterr().out)["bytes"] == 27
    with pytest.raises(SystemExit) as error:
        run(["verify", str(tmp_path / "missing")])
    assert error.value.code == 2
    assert "Cannot read bundle manifest" in capsys.readouterr().err
