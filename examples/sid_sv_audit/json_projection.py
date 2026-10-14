"""Read selected evidence from complete JSON without retaining other fields."""

import gzip
import json
from pathlib import Path

from isovar._audit_json_reader import json_events


class _ChildDecoder(json.JSONDecoder):
    """Keep the root explicit while decoding buffered child containers in C."""

    def __init__(self):
        super().__init__()
        self.first = True

    def raw_decode(self, document, index=0):
        if self.first:
            self.first = False
            raise json.JSONDecodeError("Read the root incrementally", document, index)
        return super().raw_decode(document, index)


def collection_items(path, field):
    """Yield one child of a root collection while validating the whole document.

    Arrays yield index/value pairs; objects yield key/value pairs. Retained
    memory is proportional to one selected child, not the complete result.
    Duplicate selected root fields are rejected rather than combined.
    """
    path = Path(path)
    with gzip.open(path, "rt", encoding="utf-8") if path.suffix == ".gz" else path.open(encoding="utf-8") as handle:
        events = iter(json_events(handle, decoder=_ChildDecoder()))

        def build(first):
            event, value = first
            if event == "scalar":
                return value
            if event == "start_map":
                result = {}
                for event, key in events:
                    if event == "end_map":
                        return result
                    result[key] = build(next(events))
            else:
                result = []
                for first in events:
                    if first[0] == "end_array":
                        return result
                    result.append(build(first))

        def discard(first):
            depth = int(first[0].startswith("start_"))
            while depth:
                event, _ = next(events)
                depth += int(event.startswith("start_")) - int(event.startswith("end_"))

        def children(first):
            event, value = first
            if event == "scalar":
                if not isinstance(value, (list, dict)):
                    raise ValueError("Expected JSON collection: " + field)
                yield from enumerate(value) if isinstance(value, list) else value.items()
            elif event == "start_map":
                for event, key in events:
                    if event == "end_map":
                        break
                    yield key, build(next(events))
            elif event == "start_array":
                index = 0
                for first in events:
                    if first[0] == "end_array":
                        break
                    yield index, build(first)
                    index += 1
            else:
                raise ValueError("Expected JSON collection: " + field)

        # Do not bulk-decode the root itself: seeing each root key also lets us
        # reject duplicate selected fields consistently for small documents.
        first = next(events)
        if first[0] == "start_map":
            found = False
            for event, key in events:
                if event == "end_map":
                    break
                first = next(events)
                if key == field:
                    if found:
                        raise ValueError("Duplicate JSON root collection: " + field)
                    found = True
                    yield from children(first)
                else:
                    discard(first)
            if not found:
                raise ValueError("Missing JSON root collection: " + field)
        else:
            raise ValueError("Expected JSON root object")
        if next(events, None) is not None:
            raise ValueError("Unexpected JSON event after root object")


def candidate_evidence(path, candidate):
    """Read precisely the observation and SAM rows used by one candidate."""
    ids = {w["observation"] for occurrence in candidate["occurrences"] for w in occurrence["witnesses"]}
    observations = {key: value for key, value in collection_items(path, "observations") if key in ids}
    if set(observations) != ids:
        raise ValueError("Missing candidate observations")
    records = {key for row in observations.values() for key in row["records"]}
    original = {key: value for key, value in collection_items(path, "original_records") if key in records}
    if set(original) != records:
        raise ValueError("Missing original candidate records")
    return dict(observations=observations, original_records=original)
