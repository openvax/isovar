"""Tests for :mod:`isovar.json_files`, the shared JSON reader/writer."""

from isovar.json_files import read_json, write_json


def test_write_json_sort_keys_and_trailing_newline(tmp_path):
    path = tmp_path / "out.json"
    write_json(path, {"b": 1, "a": 2}, sort_keys=True)
    assert path.read_text() == '{\n  "a": 2,\n  "b": 1\n}\n'
    write_json(path, {"b": 1, "a": 2})
    assert path.read_text().startswith('{\n  "b"')
    assert read_json(path) == {"b": 1, "a": 2}
