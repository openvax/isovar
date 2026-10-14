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
