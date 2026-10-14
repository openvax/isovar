"""Streaming must preserve original bytes, values and publication guarantees."""

from copy import deepcopy
import gzip
from hashlib import sha256
import itertools
import json
import tracemalloc

import pytest

from isovar._audit_json import canonical_chunks, comparison_identity, identity, normalize_json
from isovar.audit_bundle import publish_json
from examples.sid_sv_audit.inventory import read_json, write_json


VALUES = [None, True, False, 0, -1, 2**100, 0.0, -0.0, 1.0, 1e-200, 1e200,
          "", "quoted\"\\\t\n\x00", "\u00e9\U0001f600\ud800", [], {},
          [1, (2, 3), {"b": False, "a": None}], {1: "a", 2: (True, False)}]


@pytest.mark.parametrize("value", VALUES)
def test_streaming_canonical_bytes_and_identity_match_released_format(value):
    expected = json.dumps(value, sort_keys=True, separators=(",", ":"), allow_nan=False).encode()
    assert b"".join(canonical_chunks(value)) == expected
    assert identity(value) == sha256(expected).hexdigest()
    normalized = normalize_json(deepcopy(value))
    assert normalized == json.loads(expected)
    assert identity(normalized) == identity(json.loads(expected))


@pytest.mark.parametrize("compressed", [False, True])
@pytest.mark.parametrize("large", [False, True])
def test_streamed_writes_are_byte_identical_and_reuse_old_checkpoints(tmp_path, compressed, large):
    value = {"rows": VALUES * (10000 if large else 1)}
    expected = json.dumps(value, sort_keys=True, separators=(",", ":"), allow_nan=False).encode() + b"\n"
    if compressed:
        expected = gzip.compress(expected, mtime=0)
    path = tmp_path / ("result.json.gz" if compressed else "result.json")
    write_json(path, value)
    assert path.read_bytes() == expected
    path.unlink()
    path.write_bytes(expected)  # A checkpoint made by the previous writer.
    write_json(path, value)
    assert path.read_bytes() == expected
    assert read_json(path) == json.loads(gzip.decompress(expected) if compressed else expected)


@pytest.mark.parametrize("compressed", [False, True])
def test_oversized_strings_between_small_values_preserve_every_byte(tmp_path, compressed):
    value = [{"before": [None, False, -0.0]}, "ACGT" * 100000,
             {"between": "quoted\"\\\n"}, "\U0001f600\ud800" * 20000, [], "after"]
    canonical = json.dumps(value, sort_keys=True, separators=(",", ":"), allow_nan=False).encode()
    assert identity(value) == sha256(canonical).hexdigest()
    expected = canonical + b"\n"
    if compressed:
        expected = gzip.compress(expected, mtime=0)
    path = tmp_path / ("result.json.gz" if compressed else "result.json")
    write_json(path, value)
    assert path.read_bytes() == expected


@pytest.mark.parametrize("bad", [float("nan"), float("inf"), object()])
@pytest.mark.parametrize("writer", [publish_json, write_json])
def test_failed_serialization_cleans_up_without_publication(tmp_path, writer, bad):
    path = tmp_path / "checkpoint.json.gz"
    with pytest.raises((TypeError, ValueError)):
        writer(path, {"late_invalid_value": ["x" * 100000, bad]})
    assert list(tmp_path.iterdir()) == []


def test_in_place_normalization_preserves_shared_evidence_and_json_key_semantics():
    shared = {"placements": [(1, 2), (3, 4)]}
    result = {"a": [shared, shared], "b": shared}
    expected = json.loads(json.dumps(result))
    assert normalize_json(result) is result
    assert result == expected
    assert result["a"][0] is result["a"][1] is result["b"] is shared
    keys = {None: (1, 2)}
    assert normalize_json(keys) == {"null": [1, 2]}
    with pytest.raises(TypeError):
        normalize_json({1: "one", "2": "two"})  # sort_keys rejected this too.


def test_comparison_fingerprint_matches_json_equality_for_numbers_and_structure():
    values = [json.loads(json.dumps(value)) for value in VALUES]
    values += [1, [True], [1.0], {"x": -0.0}, {"x": 0},
               {"x": 1, "y": [False]}, {"y": [0.0], "x": True},
               2**53, float(2**53), 2**53 + 1, "null", [[], []], [[]],
               {"ab": "c"}, {"a": "bc"}]
    fingerprints = [comparison_identity(value) for value in values]
    for i, j in itertools.product(range(len(values)), repeat=2):
        assert (fingerprints[i] == fingerprints[j]) == (values[i] == values[j])


def test_hashing_and_writing_do_not_allocate_an_entire_serialized_document(tmp_path):
    # Existing evidence is excluded from the measurement: observe only the
    # operation's additional memory. Its canonical JSON is over 10 MiB.
    rows = [{"sam": "ACGT" * 256, "record": i} for i in range(12000)]
    for operation in (lambda: identity(rows), lambda: write_json(tmp_path / "large.json.gz", rows)):
        tracemalloc.start()
        try:
            operation()
            _, peak = tracemalloc.get_traced_memory()
        finally:
            tracemalloc.stop()
        assert peak < 2 * 1024 * 1024
