import gzip
import json
import tracemalloc

import pytest

from examples.sid_sv_audit.json_projection import candidate_evidence, collection_items


@pytest.mark.parametrize("compressed", [False, True])
@pytest.mark.parametrize("chunk_size", [1, 7, 65536])
def test_selected_collections_match_full_document(tmp_path, monkeypatch, compressed, chunk_size):
    from isovar import _audit_json_reader
    monkeypatch.setattr(_audit_json_reader, "_CHUNK_SIZE", chunk_size)
    path = tmp_path / ("result.json.gz" if compressed else "result.json")
    value = dict(before=[{"unused": "text"}], paths=[{"dna": "ATG", "nested": [1, True, None]}, {"dna": "é\\\""}],
                 observations={"one": {"records": ["r"]}, "two": {}}, after={"unused": []})
    encoded = json.dumps(value).encode()
    path.write_bytes(gzip.compress(encoded) if compressed else encoded)
    assert list(collection_items(path, "paths")) == list(enumerate(value["paths"]))
    assert dict(collection_items(path, "observations")) == value["observations"]


@pytest.mark.parametrize("text", [
    '{"paths":[],"ignored":[1,]}', '{"paths":[{}]', '{"paths":[]} false',
    '{"paths":[],"paths":[]}', '{"other":[]}', '{"paths":1}', '[1,2]',
])
def test_invalid_complete_document_is_rejected(tmp_path, text):
    path = tmp_path / "result.json"
    path.write_text(text)
    with pytest.raises(ValueError):
        list(collection_items(path, "paths"))


def test_candidate_projection_works_independently_of_root_field_order(tmp_path):
    path = tmp_path / "result.json"
    value = dict(original_records={"r1": "original sam", "r2": "background sam"},
                 paths=[{"unused": [1, 2]}], observations={"o1": {"records": ["r1"]}, "o2": {"records": ["r2"]}})
    path.write_text(json.dumps(value))
    candidate = dict(occurrences=[dict(witnesses=[dict(observation="o1")])])
    assert candidate_evidence(path, candidate) == dict(observations={"o1": value["observations"]["o1"]},
                                                     original_records={"r1": "original sam"})
    del value["original_records"]["r1"]
    path.write_text(json.dumps(value))
    with pytest.raises(ValueError, match="Missing original candidate records"):
        candidate_evidence(path, candidate)
    del value["observations"]["o1"]
    path.write_text(json.dumps(value))
    with pytest.raises(ValueError, match="Missing candidate observations"):
        candidate_evidence(path, candidate)


def test_background_collection_is_not_retained(tmp_path):
    path = tmp_path / "result.json"
    with path.open("w") as handle:
        handle.write('{"unused":[')
        for i in range(30000):
            handle.write(("," if i else "") + json.dumps({"sequence": "ACGT" * 25, "index": i}))
        handle.write('],"paths":[{"sequence":"ATG"}]}')
    tracemalloc.start()
    try:
        assert list(collection_items(path, "paths")) == [(0, {"sequence": "ATG"})]
        _, peak = tracemalloc.get_traced_memory()
    finally:
        tracemalloc.stop()
    assert peak < 2 * 1024**2
