"""Read selected evidence from complete JSON."""

import gzip
import json
from pathlib import Path

from isovar.sv_rna_orf_export import WitnessPool


def read_evidence_json(handle):
    """Decode complete evidence, sharing equal nested witness values only."""
    keys, pool = {}, WitnessPool()
    witness_fields = {"observation", "observation_interval", "original_query_interval", "junction_links"}

    def normalize(row):
        row = {(keys.setdefault(k, k) if len(keys) < 1024 else keys.get(k, k)): v
               for k, v in row.items()}
        if witness_fields <= row.keys():
            return {k: pool.share(v) for k, v in row.items()}
        return row

    result = json.load(handle, object_hook=normalize)
    pool.strings.clear()
    pool.containers.clear()
    return result


def collection_items(path, field):
    """Yield one child of a root collection while validating the whole document.

    Arrays yield index/value pairs; objects yield key/value pairs.
    Duplicate selected root fields are rejected rather than combined.
    """
    path = Path(path)
    last_pairs = [None]

    def track_pairs(pairs):
        last_pairs[0] = pairs
        return dict(pairs)

    with gzip.open(path, "rt", encoding="utf-8") if path.suffix == ".gz" else path.open(encoding="utf-8") as handle:
        document = json.load(handle, object_pairs_hook=track_pairs)
    if not isinstance(document, dict):
        raise ValueError("Expected JSON root object")
    if sum(1 for key, _ in (last_pairs[0] or ()) if key == field) > 1:
        raise ValueError("Duplicate JSON root collection: " + field)
    if field not in document:
        raise ValueError("Missing JSON root collection: " + field)
    value = document[field]
    if not isinstance(value, (list, dict)):
        raise ValueError("Expected JSON collection: " + field)
    yield from enumerate(value) if isinstance(value, list) else value.items()


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
