"""Tests for :mod:`isovar.json_files`, the shared JSON reader/writer."""

import gzip
import json
import pytest

from isovar.json_files import read_json, write_json


def test_write_json_sort_keys_and_trailing_newline(tmp_path):
    path = tmp_path / "out.json"
    write_json(path, {"b": 1, "a": 2}, sort_keys=True)
    assert path.read_text() == '{\n  "a": 2,\n  "b": 1\n}\n'
    write_json(path, {"b": 1, "a": 2})
    assert path.read_text().startswith('{\n  "b"')
    assert read_json(path) == {"b": 1, "a": 2}


def test_write_json_canonical_plain_and_gzip_and_streaming_reader(tmp_path):
    plain = tmp_path / "canonical.json"
    compressed = tmp_path / "canonical.json.gz"
    payload = {"b": [1, 2], "a": "value"}
    write_json(plain, payload, canonical=True)
    write_json(compressed, payload, canonical=True)
    assert plain.read_text() == '{"a":"value","b":[1,2]}\n'
    assert gzip.decompress(compressed.read_bytes()) == b'{"a":"value","b":[1,2]}\n'
    assert read_json(plain, loader=json.load) == payload
    assert read_json(compressed, loader=json.load) == payload
    with pytest.raises(ValueError):
        write_json(plain, {"bad": float("nan")}, canonical=True)
    assert read_json(plain) == payload
    assert sorted(p.name for p in tmp_path.iterdir()) == ["canonical.json", "canonical.json.gz"]
