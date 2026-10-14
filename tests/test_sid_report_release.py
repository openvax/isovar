import gzip
import json
from pathlib import Path

import pytest

from examples import prepare_sid_report_release as release
from examples.sid_sv_audit.inventory import digest, read_json, write_json


def report(tmp_path):
    source = tmp_path / "source"
    (source / "events").mkdir(parents=True)
    write_json(source / "events.json.gz", dict(candidates=[dict(witnesses=["original-read", "other-read"])]))
    (source / "README.md").write_text("[Complete ledger](events.json.gz)\n")
    (source / "events" / "event.md").write_text("[Complete ledger](../events.json.gz)\n[Other](other.json.gz)\n")
    (source / "support.tsv").write_text("reads\n2\n")
    pins = {str(p.relative_to(source)): digest(p) for p in source.rglob("*") if p.is_file()}
    write_json(source / "files.json", pins)
    return source


def test_complete_release_asset_and_distribution_resolve_every_checksum(tmp_path):
    source = report(tmp_path)
    before = {str(p.relative_to(source)): p.read_bytes() for p in source.rglob("*") if p.is_file()}
    destination, assets = tmp_path / "distributed", tmp_path / "assets"
    download = release.prepare(source, destination, assets, "1.54.0")
    asset = assets / download["filename"]
    assert asset.read_bytes() == before["events.json.gz"]
    assert json.loads(gzip.decompress(asset.read_bytes()))["candidates"][0]["witnesses"] == ["original-read", "other-read"]
    assert not (destination / "events.json.gz").exists()
    assert download["url"] in (destination / "README.md").read_text()
    assert download["url"] in (destination / "events/event.md").read_text()
    assert "[Other](other.json.gz)" in (destination / "events/event.md").read_text()
    inventory = read_json(destination / "downloads.json")
    assert inventory["files"] == {"events.json.gz": download}
    assert inventory["source_files_sha256"] == digest(source / "files.json")
    for name, checksum in read_json(destination / "files.json").items():
        assert digest(asset if name == "events.json.gz" else destination / name) == checksum
    assert before == {str(p.relative_to(source)): p.read_bytes() for p in source.rglob("*") if p.is_file()}
    # An identical pre-existing asset is reusable; the source stays immutable.
    assert release.prepare(source, tmp_path / "second", assets, "1.54.0") == download
    with pytest.raises(ValueError, match="already exists"):
        release.prepare(source, destination, assets, "1.54.0")


@pytest.mark.parametrize("change", ["checksum", "extra", "asset", "size", "version"])
def test_invalid_release_preparation_never_publishes_a_report(tmp_path, monkeypatch, change):
    source = report(tmp_path)
    destination, assets = tmp_path / "distributed", tmp_path / "assets"
    version = "1.54.0"
    if change == "checksum":
        (source / "support.tsv").write_text("changed")
    elif change == "extra":
        (source / "unexpected").write_text("untracked")
    elif change == "asset":
        assets.mkdir()
        (assets / "sid-neoorf-events-1.54.0.json.gz").write_bytes(b"different asset")
    elif change == "size":
        monkeypatch.setattr(release, "_MAX_ASSET_BYTES", (source / "events.json.gz").stat().st_size)
    else:
        version = "../../elsewhere"
    with pytest.raises(ValueError):
        release.prepare(source, destination, assets, version)
    assert not destination.exists()
    assert (source / "events.json.gz").exists()
    if change == "asset":
        assert (assets / "sid-neoorf-events-1.54.0.json.gz").read_bytes() == b"different asset"


@pytest.mark.parametrize("layout", ["source_output", "source_assets", "same_output", "nested_assets"])
def test_overlapping_locations_are_rejected_before_mutation(tmp_path, layout):
    source = report(tmp_path)
    before = {str(p.relative_to(source)): p.read_bytes() for p in source.rglob("*") if p.is_file()}
    output, assets = tmp_path / "output", tmp_path / "assets"
    if layout == "source_output":
        output = source / "output"
    elif layout == "source_assets":
        assets = source / "assets"
    elif layout == "same_output":
        assets = output
    else:
        assets = output / "assets"
    with pytest.raises(ValueError, match="separate locations"):
        release.prepare(source, output, assets, "1.54.0")
    assert not output.exists() and not assets.exists()
    assert before == {str(p.relative_to(source)): p.read_bytes() for p in source.rglob("*") if p.is_file()}


@pytest.mark.parametrize("change", ["file", "manifest_bytes", "new_file"])
def test_source_drift_during_copy_prevents_distribution_publication(tmp_path, monkeypatch, change):
    source = report(tmp_path)
    copy = release.shutil.copyfile

    def drift(src, dest):
        result = copy(src, dest)
        if Path(src).name == "events.json.gz":
            if change == "file":
                (source / "support.tsv").write_text("changed support")
            elif change == "manifest_bytes":
                manifest = source / "files.json"
                manifest.write_text(manifest.read_text() + " ")
            else:
                (source / "new.txt").write_text("unexpected")
        return result

    monkeypatch.setattr(release.shutil, "copyfile", drift)
    with pytest.raises(ValueError, match="changed during preparation"):
        release.prepare(source, tmp_path / "distributed", tmp_path / "assets", "1.54.0")
    assert not (tmp_path / "distributed").exists()
    assert not list(tmp_path.glob(".report-release-*"))
    assert not list((tmp_path / "assets").glob(".ledger-*"))
