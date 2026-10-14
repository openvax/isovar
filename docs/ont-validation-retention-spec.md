# Bound complete-result validation (#500)

Both frozen ONT reconstructions completed. The projected forward verifier checked
all 200 path fingerprints and 699 native ORF translations, then exceeded its
24 GiB sampled physical-footprint guard while retaining subsequent evidence.
Preserve that attempt, its code and logs; do not publish either successor yet.

Restore exact sharing of nested witness values (observation IDs, interval
arrays and junction links) by their complete canonical JSON values. Whole rows
can differ by path-relative coordinates and must remain distinct. Preserve every occurrence and observation identity,
including different intervals, links, numeric types and signed zero. Drop full
path temporaries before decoding the next path and release the interning index
before observation/SAM loading. Keep validating every complete path and every
original SAM record; keep the producing files and all scientific limits frozen.

The standard ORF exporter also round-trips entire occurrence objects through
JSON, recreating shared witness dictionaries. Share equal nested values inside the exported witness copies, keeping them
independent from the input. Deduplicate and order occurrence values by streamed
canonical metadata and individual witness values. Verify complete values even
on a checksum match, preserving the old full canonical-JSON order, including
array-prefix ordering. Normalize metadata to the
same JSON values without materializing a complete occurrence string.

Build the unchanged cell/lineage lookup maps once per export. Repeated uses of
one observation key need one metadata lookup: every occurrence remains in the
output, while the existing unique read/fragment/support identities and flags
remain identical. This changes representation and repeated work, not support
units, original read grouping or scientific interpretation.

Compare complete old/new exports and summaries on existing small and short-read
fixtures; test canonical ordering, duplicates, alias isolation from the input,
and distinct witness values. Record the updated validator/exporter identity
separately from the frozen producer, then rerun both full independent checks
under the same 4-hour, 16 GiB RSS, 24 GiB physical-footprint and 8 GiB disk limits.
Rebuild the separately pinned installed replay wheel and validate its full source
identity and offline fixture before continuing staged archival and full replay.

Summary and report readers must retain the same sharing when reopening saved
evidence. Decode incrementally and intern equal nested values only in witness
rows, including rows decoded inside buffered containers. Keep arbitrary JSON
values and every witness row unchanged; test both buffered and streamed paths
against the standard decoder. The frozen producing helpers remain untouched.

The shared forward verifier completed under its limits. Reverse completed every
scientific check and wrote its receipts, then exceeded the sampled RSS limit
during process teardown. Preserve that failed process record and its outputs.
Before loading original records, retain cell/lineage rows only for the original
identities used by selected candidate observations. Release projected paths and
metadata after exporting the complete summary. Re-run reverse with the same
guards and require its summary bytes to match the preceding computed summary.
Retain the successful forward proof through an explicit per-view verification
inventory; do not rerun or relabel it as produced by the new verifier.

Before any publication stage, sample the entire isolated stage process group
(#502). Report creation launches subprocesses; sampling only the coordinator
misses their memory. Sum live RSS and physical-footprint measurements across
that group, fail closed on a live-process measurement error, and terminate the
whole group at the unchanged limits. Test a child allocation that exceeds a
small fixture guard while its parent remains below it. Preserve the stopped
parent-only queue; no publication stage ran under that declaration.

## Resume reverse restore after reclaiming approved snapshot space

The full forward installed replay and CLI resume passed. Preserve those proofs,
both completed archive identities, and the reverse restore's failed disk-space
preflight. Resume at reverse restore after macOS reclaims user-approved local
Time Machine snapshot space; retain the same memory, time and free-space guards.
Do not repeat the successful reconstruction, replay, import or archive stages.

Record the resumed stages separately, then regenerate and archive the complete
reports. Include both preceding attempts and the continuation's declarations in
the final provenance. Update release dependency pins only after the frozen
recovery stages finish, preserving their original dependency catalogue.
