"""Packaging guardrails for isovar wheels (#678)."""

from __future__ import annotations

import builtins
import email
import zipfile
from pathlib import Path

import pytest
from setuptools import build_meta

from isovar import sid_data

REPO_ROOT = Path(__file__).resolve().parents[1]
MAX_WHEEL_BYTES = 2 * 1024 * 1024


def test_wheel_excludes_sid_reads_and_keeps_osteosarc_optional(tmp_path: Path, monkeypatch):
    dist_dir = tmp_path / "dist"
    dist_dir.mkdir()
    monkeypatch.chdir(REPO_ROOT)
    wheel_name = build_meta.build_wheel(str(dist_dir))
    wheel_path = dist_dir / wheel_name

    assert wheel_path.is_file()
    assert wheel_path.stat().st_size < MAX_WHEEL_BYTES, (
        f"Wheel {wheel_name} is {wheel_path.stat().st_size} bytes "
        f"(exceeds {MAX_WHEEL_BYTES}-byte budget)"
    )

    with zipfile.ZipFile(wheel_path) as zf:
        names = zf.namelist()
        assert not any("data/sid-reads" in name for name in names), (
            "Sid test reads must not be packaged inside the isovar wheel"
        )
        metadata_names = [name for name in names if name.endswith(".dist-info/METADATA")]
        assert len(metadata_names) == 1
        meta = email.message_from_bytes(zf.read(metadata_names[0]))

    requires = meta.get_all("Requires-Dist") or []
    osteosarc_reqs = [req for req in requires if req.lower().startswith("osteosarc")]
    assert osteosarc_reqs, "expected osteosarc in the optional [test] extra"
    assert all("extra ==" in req and "test" in req for req in osteosarc_reqs), (
        f"osteosarc must only be an optional [test] dependency, got: {osteosarc_reqs}"
    )


def test_sid_data_gives_clear_error_when_osteosarc_missing(monkeypatch):
    real_import = builtins.__import__

    def blocked_import(name, *args, **kwargs):
        if name == "osteosarc" or name.startswith("osteosarc."):
            raise ImportError("No module named 'osteosarc'")
        return real_import(name, *args, **kwargs)

    monkeypatch.setattr(builtins, "__import__", blocked_import)
    with pytest.raises(ImportError, match=r"isovar\[test\]"):
        sid_data.bundle()
