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
@pytest.mark.parametrize("chunk_size", [1, 7, 65536])
def test_evidence_reader_shares_nested_witness_values_on_bulk_and_streamed_paths(monkeypatch, chunk_size):
    import io
    from examples.sid_sv_audit.json_projection import read_evidence_json
    from isovar import _audit_json_reader
    monkeypatch.setattr(_audit_json_reader, "_CHUNK_SIZE", chunk_size)

    witness = dict(observation="same-read", observation_interval=[1, 19],
                   original_query_interval=[31, 49], junction_links=[dict(query_interval=[7, 7])])
    for padding in ("", "x" * 70000 if chunk_size == 65536 else "x" * 30):
        value = dict(candidates=[dict(padding=padding, witnesses=[witness, witness]),
                                 dict(witnesses=[dict(witness, original_query_interval=[32, 50])])],
                     other=[False, 0, 0.0, -0.0, None, "é"])
        text = json.dumps(value, sort_keys=True)
        actual = read_evidence_json(io.StringIO(text))
        assert json.dumps(actual, sort_keys=True) == text
        first, second = actual["candidates"][0]["witnesses"]
        third = actual["candidates"][1]["witnesses"][0]
        assert first is not second
        assert first["observation_interval"] is second["observation_interval"] is third["observation_interval"]
        assert first["junction_links"] is third["junction_links"]
        assert first["original_query_interval"] is not third["original_query_interval"]


def test_evidence_reader_validates_complete_json_and_preserves_scalar_roots():
    import io
    from examples.sid_sv_audit.json_projection import read_evidence_json

    for text in ("null", "false", "12", '"string"', "[]", "{}", '[1, {"a": 2, "a": 3}]'):
        assert read_evidence_json(io.StringIO(text)) == json.loads(text)
    for text in ("", '{"a":', "[1,]", "{} {}", '{"a": 1,}'):
        with pytest.raises(ValueError):
            read_evidence_json(io.StringIO(text))
