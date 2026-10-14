# Stream CLI checkpoint verification (#498)

The audit CLI returns a status and checksum, so resumed CLI runs should not
materialize the full reconstruction. Preserve the public replay API's existing
full-result return contract, the checkpoint JSON values, and all validation gates.

For generated checkpoints, consume the result's JSON events once to compute
its canonical SHA-256 and complete semantic comparison fingerprint. Retain
only request/status/checksum metadata. Validate the entire document, including
trailing content, before returning success. Canonical object key order permits
streamed hashing; fall back to the existing materialized validator for unsorted
large objects or duplicate root keys, preserving their current JSON semantics.
Small buffered objects remain bounded and can be canonicalized directly.

Keep expected-result comparison streamed, verify bundle/software/request pins,
and reject changed result checksums or comparison flags exactly as before.
New CLI replays still reconstruct normally, then release the full result before
the next run. New checkpoints use deterministic gzip JSON to avoid allocating
an additional uncompressed reconstruction on disk. Existing plain checkpoints
remain readable and byte-identical on resume. Reject ambiguous plain/gzip pairs
rather than picking one silently. Add no command-line switches or new public
tuning parameters.

Test full API/CLI equivalence, unchanged checkpoint bytes, tampered request and
result data, comparison changes, malformed trailing JSON, noncanonical objects,
duplicates, numeric equality and large arrays. Measure summary-only resume
memory and verify the installed wheel offline before release.
