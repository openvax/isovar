"""Shared fixture bytes, native identities and explicit offline acquisition."""

from collections import Counter
from copy import deepcopy
from hashlib import sha256
from io import BytesIO
import json
from pathlib import Path
import shutil
import subprocess
import sys

import pysam
import pytest
import requests

from isovar.osteosarc_data import DEFAULT_MANIFEST, acquire_dataset, load_manifest, main, verify_dataset


CORPUS = Path(__file__).parent / "data/osteosarc/expansion/corpus"


@pytest.fixture
def osteosarc():
    if sys.version_info < (3, 10):
        pytest.skip("Optional acquisition uses osteosarc, which requires Python 3.10+")
    return pytest.importorskip("osteosarc")


@pytest.fixture
def no_network(monkeypatch):
    monkeypatch.setattr(subprocess, "run", lambda *a, **kw: pytest.fail("Unexpected network/acquisition subprocess"))
    monkeypatch.setattr(requests.Session, "send", lambda *a, **kw: pytest.fail("Unexpected HTTP request"))


@pytest.fixture
def small_manifest(tmp_path):
    manifest = deepcopy(load_manifest())
    manifest["cases"] = manifest["cases"][:1]
    manifest["assets"] = manifest["assets"][:2]
    path = tmp_path / "manifest.json"
    path.write_text(json.dumps(manifest))
    return path


def test_shared_manifest_preserves_original_native_cases_and_asset_bytes():
    manifest = load_manifest()
    original = json.loads((CORPUS / "manifest.json").read_text())
    assert manifest["data_version"] == "vaccine-rna-v1"
    assert len(manifest["cases"]) == len(original["cases"]) == 49
    assert len(manifest["assets"]) == 98
    assert sum(a["size_bytes"] for a in manifest["assets"]) == 4126730
    assert sha256((CORPUS / "manifest.json").read_bytes()).hexdigest() == manifest["upstream"]["manifest_sha256"]
    for case, source in zip(manifest["cases"], original["cases"]):
        for key in ("case_id", "variant", "selection", "source_id", "source_url",
                    "source_region_sha256", "source_receipt_sha256", "selection_row_sha256"):
            assert case[key] == source[key]
        assert case["selected_record_count"] == len(source["selected_records"])
    assert {c["variant"]["assembly"] for c in manifest["cases"]} == {"GRCh37", "GRCh38"}
    for asset in manifest["assets"]:
        raw = (CORPUS / asset["filename"]).read_bytes()
        assert sha256(raw).hexdigest() == asset["sha256"] and len(raw) == asset["size_bytes"]


def test_offline_import_and_export_preserve_every_original_record(tmp_path, osteosarc, no_network):
    cache = tmp_path / "cache"
    paths = acquire_dataset(cache_root=cache, import_corpus=CORPUS, offline=True)
    identities = {name: (p.stat().st_ino, p.stat().st_mtime_ns) for name, p in paths.items()}
    output = acquire_dataset(tmp_path / "export", cache_root=cache, offline=True)
    assert verify_dataset(output) == output
    for case in load_manifest()["cases"]:
        with pysam.AlignmentFile(str(CORPUS / case["bam"])) as original, \
                pysam.AlignmentFile(str(output / case["bam"])) as exported:
            assert str(exported.header) == str(original.header)
            assert Counter(r.to_string() for r in original) == Counter(r.to_string() for r in exported)
    assert {name: (p.stat().st_ino, p.stat().st_mtime_ns) for name, p in paths.items()} == identities
    # An identical existing export needs neither the cache nor osteosarc acquisition.
    assert acquire_dataset(output, cache_root=tmp_path / "absent", offline=True) == output
    assert not (tmp_path / "absent").exists()


def test_reuses_vaxrank_datacache_objects_without_url_receipts(tmp_path, osteosarc, small_manifest, no_network):
    from datacache import Cache
    shared = tmp_path / "shared"
    cache = Cache("openvax", cache_root=shared / "objects/sha256")
    objects = []
    for asset in load_manifest(small_manifest)["assets"]:
        path = Path(cache.local_path(filename=asset["sha256"] + "".join(Path(asset["filename"]).suffixes)))
        path.parent.mkdir(parents=True, exist_ok=True)
        shutil.copyfile(CORPUS / asset["filename"], path)
        objects.append((path, path.stat().st_ino, path.stat().st_mtime_ns))
    assert not (shared / "osteosarc").exists()
    acquire_dataset(tmp_path / "export", manifest_path=small_manifest, cache_root=shared, offline=True)
    assert all((p.stat().st_ino, p.stat().st_mtime_ns) == (inode, modified) for p, inode, modified in objects)


def test_missing_cache_fails_offline_before_export(tmp_path, osteosarc, small_manifest, no_network):
    with pytest.raises(ValueError, match="Not cached"):
        acquire_dataset(tmp_path / "export", manifest_path=small_manifest, cache_root=tmp_path / "cache", offline=True)
    assert not (tmp_path / "export").exists()


def test_corrupt_cache_and_modified_exports_are_not_silently_repaired(tmp_path, osteosarc, small_manifest, no_network):
    cache = tmp_path / "cache"
    output = acquire_dataset(tmp_path / "export", manifest_path=small_manifest, cache_root=cache,
                             import_corpus=CORPUS, offline=True)
    asset = load_manifest(small_manifest)["assets"][0]
    cached = cache / "objects/sha256" / (asset["sha256"] + ".bam")
    cached.write_bytes(b"x" * cached.stat().st_size)
    with pytest.raises(ValueError, match="checksum mismatch"):
        acquire_dataset(manifest_path=small_manifest, cache_root=cache, offline=True)
    (output / asset["filename"]).write_bytes(b"modified")
    with pytest.raises(ValueError, match="checksum mismatch"):
        acquire_dataset(output, manifest_path=small_manifest, cache_root=cache, repair_cache=True)


@pytest.mark.parametrize("corrupt_download", [False, True])
@pytest.mark.parametrize("corruption", ["object", "url_receipt"])
def test_explicit_repair_checks_downloaded_bytes(tmp_path, osteosarc, small_manifest, monkeypatch,
                                               corrupt_download, corruption, no_network):
    cache = tmp_path / "cache"
    paths = acquire_dataset(manifest_path=small_manifest, cache_root=cache, import_corpus=CORPUS, offline=True)
    asset = load_manifest(small_manifest)["assets"][0]
    if corruption == "object":
        paths[asset["filename"]].write_bytes(b"corrupt")
    else:
        paths[asset["filename"]].unlink()
        wrong = tmp_path / "wrong.bam"
        wrong.write_bytes(b"x" * asset["size_bytes"])
        osteosarc.Cache(cache).import_file(wrong, asset["url"])
    calls = []

    def download(session, request, **kwargs):
        assert request.url == asset["url"]
        assert request.method in ("HEAD", "GET")
        calls.append(request.method)
        response = requests.Response()
        response.status_code = 200
        response.url, response.request = request.url, request
        response.headers.update({"Content-Length": str(asset["size_bytes"]), "ETag": '"fixture"'})
        payload = b"x" * asset["size_bytes"] if corrupt_download else (CORPUS / asset["filename"]).read_bytes()
        response.raw = BytesIO(payload if request.method == "GET" else b"")
        return response

    # Exercise the real osteosarc/datacache download and validation path;
    # replace only the HTTP transport, keeping the test independent of curl.
    monkeypatch.setattr(requests.Session, "send", download)
    if corrupt_download:
        with pytest.raises(ValueError, match="(?i)sha-?256 mismatch"):
            acquire_dataset(tmp_path / "export", manifest_path=small_manifest, cache_root=cache, repair_cache=True)
        assert not (tmp_path / "export").exists()
    else:
        output = acquire_dataset(tmp_path / "export", manifest_path=small_manifest, cache_root=cache, repair_cache=True)
        assert verify_dataset(output, small_manifest) == output
    assert calls.count("GET") == 1
    assert calls.count("HEAD") >= 1


def test_cli_offline_failure_is_actionable(tmp_path, osteosarc, small_manifest, no_network, capsys):
    with pytest.raises(SystemExit) as error:
        main(["--manifest", str(small_manifest), "--cache-root", str(tmp_path / "cache"), "--offline"])
    assert error.value.code == 1
    assert "Osteosarc fixtures unavailable: Not cached" in capsys.readouterr().err


@pytest.mark.parametrize("filename", ["../outside.bam", "manifest.json", "/outside.bam"])
def test_manifest_rejects_unsafe_asset_paths(tmp_path, filename):
    manifest = load_manifest()
    manifest["assets"][0]["filename"] = filename
    path = tmp_path / "bad.json"
    path.write_text(json.dumps(manifest))
    with pytest.raises(ValueError, match="basenames"):
        load_manifest(path)
