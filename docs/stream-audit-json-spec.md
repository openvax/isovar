# Stream audit JSON (#488)

Keep canonical JSON bytes, SHA-256 identities, complete evidence, atomic
publication and replay comparison semantics unchanged while avoiding complete
serialized copies. Share canonical encoding/hash helpers between the installed
audit API and the Sid runner. Stream plain and deterministic gzip writes into
temporary files, compare existing checkpoints in chunks, and publish only after
serialization and fsync succeed. Read gzip through a file stream so compressed
input and a separate decompressed byte buffer are not retained together.

Replay fingerprints the expected parsed result and releases it before running
the engine. Normalize native tuples to JSON lists in place while retaining
shared containers. Compare all JSON values with a deterministic fingerprint
that preserves existing Python JSON equality (including numeric equality),
separate from the unchanged canonical checksum used for checkpoint integrity.
Reuse and mismatch diagnostics keep their current behavior.

Sid import validates each document's pins before releasing its parsed payload.
Retain acquisition receipts and checkpoint provenance, not complete reports or
candidate arrays. Build each replay plan from one reconstruction at a time;
the original files remain copied verbatim and every declared pin is validated.

The standard JSON decoder still needs the decoded text and one parsed object;
this change does not claim bounded-memory parsing. Measure that independently
from hashing, writing and replay. Check canonical/gzip byte equivalence, old
checkpoint reuse, interruption/racing writers, numeric/container edge cases,
whole-result replay and peak allocation on representative evidence.

Ship a versioned PR after lint, full tests, four-version CI and wheel replay.
Then carry the shipped change into the existing ONT recovery PR, preserve prior
attempts, rerun both orientations with declared budgets and current input pins,
and finish independent checks, reports, backup/restore and installed replay.
