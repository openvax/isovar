import gc
import io
import json
import tracemalloc

import pytest

from isovar import _audit_json_reader
from isovar._audit_checkpoint import NonCanonicalCheckpoint, checkpoint_summary
from isovar._audit_json import canonical, comparison_identity, identity
from isovar.audit_bundle import publish_json


@pytest.mark.parametrize("chunk_size", [1, 17, 65536])
@pytest.mark.parametrize("result", [
    None, [], {}, [True, 1, 1.0, -0.0, 2**100, 0.125],
    {"a": ["é", "\\\"", "\ud800", {"x": []}], "z": {"b": False}},
])
def test_streamed_checkpoint_identities_match_full_values(monkeypatch, chunk_size, result):
    monkeypatch.setattr(_audit_json_reader, "_CHUNK_SIZE", chunk_size)
    checkpoint = dict(request={"run_id": "example"}, result=result,
                      result_sha256=identity(result), comparison="matched")
    metadata, checksum, fingerprint = checkpoint_summary(io.StringIO(canonical(checkpoint).decode()))
    assert metadata == {k: v for k, v in checkpoint.items() if k != "result"}
    assert checksum == identity(result)
    assert fingerprint == comparison_identity(result)


@pytest.mark.parametrize("text", [
    '{"result":{"z":1,"a":2}}', '{"result":{"a":1,"a":2}}',
    '{"result":[],"result":{}}',
])
def test_noncanonical_large_objects_request_full_reader(monkeypatch, text):
    monkeypatch.setattr(_audit_json_reader, "_CHUNK_SIZE", 1)
    with pytest.raises(NonCanonicalCheckpoint):
        checkpoint_summary(io.StringIO(text))


@pytest.mark.parametrize("text", [
    '{"result":[]} false', '{"result":[1,]}', '{"result":[]', '[]', '',
])
def test_streamed_checkpoint_validates_entire_document(text):
    with pytest.raises(ValueError):
        checkpoint_summary(io.StringIO(text))


def test_checkpoint_missing_result_is_distinct_from_literal_null():
    assert checkpoint_summary(io.StringIO('{"comparison":"matched"}')) == ({"comparison": "matched"}, None, None)
    assert checkpoint_summary(io.StringIO('{"result":null}'))[1:] == (identity(None), comparison_identity(None))


def test_noncanonical_buffered_objects_preserve_standard_json_values():
    text = '{"result":{"z":1,"a":2,"a":3},"request":{"old":0},"request":{"new":1}}'
    value = json.loads(text)
    metadata, checksum, comparison = checkpoint_summary(io.StringIO(text))
    assert metadata == {"request": {"new": 1}}
    assert checksum == identity(value["result"])
    assert comparison == comparison_identity(value["result"])


def test_streamed_large_array_does_not_retain_result_tree(tmp_path):
    result = [{"sequence": "ACGT" * 25, "i": i} for i in range(30000)]
    expected_checksum, expected_comparison = identity(result), comparison_identity(result)
    path = tmp_path / "checkpoint.json"
    publish_json(path, dict(result=result, comparison="matched", result_sha256=expected_checksum, request={}))
    del result
    gc.collect()
    tracemalloc.start()
    try:
        with path.open() as handle:
            metadata, checksum, comparison = checkpoint_summary(handle)
        _, peak = tracemalloc.get_traced_memory()
    finally:
        tracemalloc.stop()
    assert checksum == expected_checksum and comparison == expected_comparison
    assert "result" not in metadata
    assert peak < 2 * 1024**2
