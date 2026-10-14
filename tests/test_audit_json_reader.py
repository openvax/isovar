import io
import json
import math
import random
import tracemalloc

import pytest

from isovar._audit_json import comparison_identity, expected_json_identity
from isovar._audit_json_reader import read_json_stream


class ChunkOnly(io.StringIO):
    def __init__(self, text, chunk):
        super().__init__(text)
        self.chunk = chunk

    def read(self, size=-1):
        assert 0 < size <= 65536, "Reader requested a whole-document buffer"
        return super().read(min(size, self.chunk))


VALUES = [None, True, False, 123, 2**200, -0.0, 0.1, 1e-200, "", "\ud800",
          "\udfff", "é 😀\n\t\\\"", [], {}, [1, {"a": [None, -0.0, True]}],
          {"z": 1, "a.b": {"item": "x"}, "é": [2, 3]}]


@pytest.mark.parametrize("chunk", [1, 2, 7, 65536])
@pytest.mark.parametrize("value", VALUES)
def test_chunked_values_and_expected_identity(value, chunk):
    text = "\n " + json.dumps(value, ensure_ascii=True) + "\n"
    actual = read_json_stream(ChunkOnly(text, chunk))
    assert actual == json.loads(text)
    if isinstance(value, float):
        assert math.copysign(1, actual) == math.copysign(1, value)
    expected = None if value is None else comparison_identity(json.loads(text))
    assert expected_json_identity(ChunkOnly(text, chunk)) == expected


@pytest.mark.parametrize("text", [
    '{"x":1,"x":{"a":[2,3]},"a":true}', '[1e+3,-0.25,1.0e-200]',
    '{"\\ud800":1,"?":2}', '{"x":"\\ud800\\udc00"}',
    '{"x":"é\u00a0😀"}', '"' + "x" * 140000 + '"'])
def test_scalar_compatibility_and_duplicate_members(text):
    expected = json.loads(text)
    for chunk in (3, 65536):
        assert read_json_stream(ChunkOnly(text, chunk)) == expected
        assert expected_json_identity(ChunkOnly(text, chunk)) == comparison_identity(expected)


@pytest.mark.parametrize("text", [
    "", " ", "[", "{", "[1", '{"a":', '{"a" 1}', '{a:1}',
    "[1,]", '{"a":1,}', "[]{}", "true false", "null 0", "01", "-01",
    "1e", "1e+", "1.", "+1", ".1", "--2", "[tru]", "[nil]", "[1 2]",
    '"unterminated', '"bad\\x"', '"\\u12"', '"\\u12xx"', '"raw\nline"',
    "\u00a0null", "[1\u00a0,2]", '{"x":1]\n'])
def test_reject_malformed_and_truncated_json(text):
    with pytest.raises(json.JSONDecodeError):
        json.loads(text)
    for chunk in (1, 7, 65536):
        for reader in (read_json_stream, expected_json_identity):
            with pytest.raises(json.JSONDecodeError):
                reader(ChunkOnly(text, chunk))


@pytest.mark.parametrize("text", ["NaN", "Infinity", "-Infinity", "1e400", "-1e400"])
def test_preserve_standard_decoder_nonfinite_values(text):
    expected = json.loads(text)
    actual = read_json_stream(ChunkOnly(text, 1))
    assert math.isnan(actual) if math.isnan(expected) else actual == expected
    with pytest.raises((ValueError, OverflowError)):
        expected_json_identity(ChunkOnly(text, 1))


def test_streamed_comparison_does_not_retain_array_values():
    text = '["' + ('x' * 24 + '\",\"') * 99999 + 'tail"]'
    stream = ChunkOnly(text, 65536)
    tracemalloc.start()
    fingerprint = expected_json_identity(stream)
    _, peak = tracemalloc.get_traced_memory()
    tracemalloc.stop()
    assert fingerprint == comparison_identity(json.loads(text))
    assert peak < 2 * 1024**2


def test_generated_json_and_truncations():
    rng = random.Random(496)

    def value(depth):
        if depth and rng.random() < 0.7:
            children = [value(depth - 1) for _ in range(rng.randrange(6))]
            return children if rng.random() < 0.5 else {str(i): v for i, v in enumerate(children)}
        return rng.choice(VALUES[:14])

    for _ in range(50):
        text = json.dumps(value(4), ensure_ascii=bool(rng.randrange(2)))
        assert read_json_stream(ChunkOnly(text, 5)) == json.loads(text)
        for length in rng.sample(range(len(text)), min(3, len(text))):
            try:
                expected = json.loads(text[:length])
            except json.JSONDecodeError:
                with pytest.raises(json.JSONDecodeError):
                    read_json_stream(ChunkOnly(text[:length], 5))
            else:
                assert read_json_stream(ChunkOnly(text[:length], 5)) == expected


@pytest.mark.parametrize("chunk", [1, 7, 20, 65536])
def test_object_hook_matches_standard_child_before_parent_order(chunk):
    text = '{"outer":[{"x":1},{"x":2}],"other":{"y":3}}'
    seen = []

    def hook(value):
        seen.append(tuple(value))
        return ("wrapped", value)

    expected = json.loads(text, object_hook=hook)
    expected_order = list(seen)
    seen.clear()
    assert read_json_stream(ChunkOnly(text, chunk), object_hook=hook) == expected
    assert seen == expected_order
