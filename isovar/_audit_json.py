"""JSON encoding and comparison without whole-document serialized copies."""

import gzip
from hashlib import sha256
import json

from ._audit_json_reader import json_events


_ENCODER = json.JSONEncoder(sort_keys=True, separators=(",", ":"), allow_nan=False)
_CHUNK_SIZE = 64 * 1024


def canonical(value):
    """Return canonical bytes for callers that explicitly need a byte string."""
    return _ENCODER.encode(value).encode()


def canonical_chunks(value):
    """Yield canonical UTF-8 in batches; retain at most one encoded JSON token."""
    buffer = bytearray()
    for token in _ENCODER.iterencode(value):
        # ensure_ascii=True makes every character exactly one UTF-8 byte.
        for start in range(0, len(token), _CHUNK_SIZE):
            buffer.extend(token[start:start + _CHUNK_SIZE].encode())
            if len(buffer) >= _CHUNK_SIZE:
                yield bytes(buffer)
                buffer.clear()
    if buffer:
        yield bytes(buffer)


def identity(value):
    result = sha256()
    for chunk in canonical_chunks(value):
        result.update(chunk)
    return result.hexdigest()


def write_json_stream(handle, value, *, compressed=False):
    """Write canonical JSON and a newline, with the existing gzip byte format."""
    if not compressed:
        for chunk in canonical_chunks(value):
            handle.write(chunk)
        handle.write(b"\n")
        return
    with gzip.GzipFile(filename="", mode="wb", fileobj=handle, compresslevel=9, mtime=0) as output:
        for chunk in canonical_chunks(value):
            output.write(chunk)
        output.write(b"\n")
    # Python 3.11/3.12's gzip.compress(mtime=0) uses zlib's OS header byte;
    # other supported Pythons use 255. Preserve the caller runtime's old bytes.
    handle.seek(9)
    handle.write(gzip.compress(b"", mtime=0)[9:10])
    handle.seek(0, 2)


def same_file_bytes(left, right):
    with left.open("rb") as a, right.open("rb") as b:
        while True:
            chunk = a.read(_CHUNK_SIZE)
            if chunk != b.read(_CHUNK_SIZE):
                return False
            if not chunk:
                return True


def normalize_json(value):
    """Normalize native engine containers in place, preserving shared evidence.

    Match the values produced by json.loads(canonical(value)), including JSON
    object key conversion. The following canonical checksum rejects cycles,
    nonfinite numbers and unsupported leaf types before publication.
    """
    seen, tuples = set(), {}

    def visit(item):
        key = id(item)
        if isinstance(item, tuple):
            if key not in tuples:
                result = []
                tuples[key] = (item, result)  # Keep the old tuple alive too.
                result.extend(visit(child) for child in item)
            return tuples[key][1]
        if isinstance(item, (dict, list)):
            if key in seen:
                return item
            seen.add(key)
            if isinstance(item, dict):
                if any(not isinstance(k, str) for k in item):
                    converted = {}
                    for k, child in sorted(item.items()):
                        if not isinstance(k, (str, int, float, bool, type(None))):
                            raise TypeError("Unsupported JSON object key: " + repr(k))
                        converted[k if isinstance(k, str) else json.dumps(k, allow_nan=False)] = child
                    item.clear()
                    item.update(converted)
                for k in item:
                    item[k] = visit(item[k])
            else:
                for i in range(len(item)):
                    item[i] = visit(item[i])
        return item

    return visit(value)


def comparison_identity(value):
    """Fingerprint complete JSON values using their existing equality semantics.

    Python compares equal numeric values equally, including True/1/1.0 and
    signed zero. Exact integer ratios preserve that behavior without rounding
    large integers. This internal comparison is separate from canonical pins.
    """
    cached = {}

    def visit(item):
        container = isinstance(item, (dict, list))
        if container and id(item) in cached:
            return cached[id(item)]
        if not container:
            return _scalar_comparison_identity(item)
        result = sha256()
        if isinstance(item, dict):
            result.update(b"object:")
            for key in sorted(item):
                result.update(visit(key))
                result.update(visit(item[key]))
        elif isinstance(item, list):
            result.update(b"array:")
            for child in item:
                result.update(visit(child))
        fingerprint = result.digest()
        if container and len(cached) < 4096:
            cached[id(item)] = fingerprint
        return fingerprint

    return visit(value).hex()


def _scalar_comparison_identity(item):
    result = sha256()
    if item is None:
        result.update(b"null")
    elif isinstance(item, str):
        result.update(b"string:")
        if len(item) <= _CHUNK_SIZE:
            result.update(_ENCODER.encode(item).encode())
        else:
            for chunk in canonical_chunks(item):
                result.update(chunk)
    elif isinstance(item, (bool, int, float)):
        numerator, denominator = item.as_integer_ratio() if isinstance(item, float) else (int(item), 1)
        result.update(("number:%s/%s" % (numerator, denominator)).encode())
    else:
        raise TypeError("Expected normalized JSON data")
    return result.digest()


def expected_json_identity(handle):
    """Compare a complete expected document without constructing its tree.

    Object frames retain only keys and child digests for ordering/duplicate
    resolution; array frames retain a running hash. Literal null means that
    comparison was not requested, matching the existing replay convention.
    """
    return comparison_event_identity(json_events(handle, decoder=json.JSONDecoder()), null_means_missing=True)


def comparison_event_identity(events, *, null_means_missing=False):
    """Reduce complete JSON value events to the existing comparison identity."""
    stack, missing = [], False
    for event, value in events:
        if event == "start_map":
            stack.append(["map", {}, None])
            continue
        if event == "start_array":
            stack.append(["array", sha256(b"array:"), None])
            continue
        if event == "map_key":
            stack[-1][2] = value
            continue
        if event.startswith("end_"):
            kind, children, _ = stack.pop()
            if kind == "map":
                result = sha256(b"object:")
                for key in sorted(children):
                    result.update(_scalar_comparison_identity(key))
                    result.update(children[key])
                fingerprint = result.digest()
            else:
                fingerprint = children.digest()
        else:
            fingerprint = (bytes.fromhex(comparison_identity(value)) if isinstance(value, (dict, list))
                           else _scalar_comparison_identity(value))
            missing = not stack and value is None
        if not stack:
            root = fingerprint
        elif stack[-1][0] == "map":
            stack[-1][1][stack[-1][2]] = fingerprint
        else:
            stack[-1][1].update(fingerprint)
    return None if missing and null_means_missing else root.hex()
