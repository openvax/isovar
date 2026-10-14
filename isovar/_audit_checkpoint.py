"""Stream verification identities from a canonical replay checkpoint."""

from hashlib import sha256
import json

from ._audit_json import _ENCODER, canonical_chunks, comparison_event_identity, comparison_identity, identity
from ._audit_json_reader import json_events


class NonCanonicalCheckpoint(ValueError):
    """Use the full JSON reader to preserve noncanonical object semantics."""


def _subtree(first, events):
    yield first
    depth = int(first[0].startswith("start_"))
    while depth:
        item = next(events)
        depth += int(item[0].startswith("start_")) - int(item[0].startswith("end_"))
        yield item


def _value(first, events):
    event, value = first
    if event == "scalar":
        return value
    if event == "start_map":
        result = {}
        for event, key in events:
            if event == "end_map":
                return result
            result[key] = _value(next(events), events)
    else:
        result = []
        for first in events:
            if first[0] == "end_array":
                return result
            result.append(_value(first, events))


def _canonical_events(events, checksum):
    stack = []

    def encode(value):
        if isinstance(value, str) and len(value) > 65536:
            for chunk in canonical_chunks(value):
                checksum.update(chunk)
        else:
            # Buffered child containers are small. Preserve the bounded leaf
            # path above while using the C encoder for these common values.
            checksum.update(_ENCODER.encode(value).encode())

    for event, value in events:
        if event == "map_key":
            frame = stack[-1]
            if frame[2] is not None and value <= frame[2]:
                raise NonCanonicalCheckpoint("Object keys are not in canonical order")
            if frame[1]:
                checksum.update(b",")
            encode(value)
            checksum.update(b":")
            frame[1], frame[2] = frame[1] + 1, value
        elif event.startswith("end_"):
            checksum.update(b"}" if event == "end_map" else b"]")
            stack.pop()
        else:
            if stack and stack[-1][0] == "array":
                if stack[-1][1]:
                    checksum.update(b",")
                stack[-1][1] += 1
            if event.startswith("start_"):
                checksum.update(b"{" if event == "start_map" else b"[")
                stack.append([event[6:], 0, None])
            else:
                encode(value)
        yield event, value


def checkpoint_summary(handle):
    """Return metadata and both complete result identities, without its tree.

    Parser buffers and small subtrees remain bounded. Object comparison frames
    retain keys/child hashes. Unsorted large objects or duplicate root members
    request the existing materialized fallback rather than changing semantics.
    """
    events = iter(json_events(handle, decoder=json.JSONDecoder()))
    first = next(events)
    if first[0] == "scalar":
        value = first[1]
        if not isinstance(value, dict):
            raise ValueError("Replay checkpoint must be an object")
        metadata = {key: item for key, item in value.items() if key != "result"}
        checksum = identity(value["result"]) if "result" in value else None
        comparison = comparison_identity(value["result"]) if "result" in value else None
    elif first[0] == "start_map":
        metadata, seen, checksum, comparison = {}, set(), None, None
        for event, key in events:
            if event == "end_map":
                break
            if key in seen:
                raise NonCanonicalCheckpoint("Duplicate checkpoint root member")
            seen.add(key)
            first = next(events)
            if key == "result":
                result = sha256()
                comparison = comparison_event_identity(_canonical_events(_subtree(first, events), result))
                checksum = result.hexdigest()
            elif key in ("request", "result_sha256", "comparison"):
                metadata[key] = _value(first, events)
            else:
                for _ in _subtree(first, events):
                    pass
    else:
        raise ValueError("Replay checkpoint must be an object")
    if next(events, None) is not None:
        raise ValueError("Unexpected JSON event after checkpoint")
    return metadata, checksum, comparison
