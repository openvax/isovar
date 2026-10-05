# Dense SV discovery and support (#437)

Add an opt-in `dense_support` API / `--dense-support` CLI mode. Scan the entire
original alignment input into a temporary SQLite spool, retaining original SAM
records and input-scoped segment identities. Group supplementary/alternative
records before selection. Search the requested breakpoint/exon/explicit regions;
prioritize placed event joins, then breakpoint clips, over ordinary context, then select deterministic
sequence/placement representatives under `max_records`. Duplicate abundance and
coordinate order must not hide a late rare event. Keep omitted discovery groups,
failed segment builds, query/path/extension limits explicit.

Run the existing reconstruction engine on the bounded discovery BAM. Stream the
spooled original groups again to find direct junction/clip observations for the
fixed candidate paths. Reuse the existing placement, exact-interval, read-filter,
cell/UMI and signal-lineage policies. Count full-path and ORF witnesses from all
available originals, including one-read cells and missing QUAL. Keep discovery
votes separate from exhaustive direct/full-interval support; unrelated background
must never count as event support. Retain original witness records and query
intervals, and preserve ambiguous placements and junction-base alternatives.

The complete-input support pass cannot make omitted hypotheses assessed, establish
independent molecules, or demonstrate translation. Ordinary intron-competition
assessment remains scoped to discovery records and must say so. Scratch disk scales
with input size; retained support witnesses and the JSON output can scale with
actual candidate support. The default bounded API and schema stay compatible.

Verify rare late events, order/cap-independent support for a fixed candidate,
original SAM/QUAL preservation, conflicting mates, alternative paths, filters,
explicit omissions, CLI/export behavior, and the complete lint/test/CI gates.

## Verification evidence

The regression suite includes the original ATP5MG--KMT2A short-read fixture,
whose pinned peptide stays `MAQFVRNLVEKTPALVNG`, and synthetic split-read cases
with alternate sequences, reversed partner orientation, missing QUAL, conflicting
mate labels, repeated records and a visible signal parent outside search regions.
Original-read interval matching and the normal ORF exporter are reused.

A local synthetic scale check (macOS, Python 3.11.12, pysam 0.24.1) placed
13 two-piece event segments after repeated 90-base donor background records.
With `max_records=2` and assembly disabled, every run selected two discovery
records and counted all 13 event fragments. Timings measure dense reconstruction;
peak RSS covers the Python process, including input creation and imports.

| Background records | Reconstruction seconds | Peak RSS (MiB) |
|---:|---:|---:|
| 1,000 | 0.035 | 120.7 |
| 100,000 | 2.066 | 155.3 |
| 1,000,000 | 37.365 | 160.7 |

These are synthetic background-scale checks, not throughput estimates for noisy
long reads or witness-heavy events. They do not replace a rerun of the Sid pilot.
The audit runner in draft PR #433 can opt into this mode when it is ready.

Testing also reproduced the default collector's SAM-payload identity collision,
filed as [#438](https://github.com/openvax/isovar/issues/438). Dense mode raises on
a conflicting original instead of selecting an arbitrary payload. Identical
repeat records still deduplicate, and interrupted/failed spools are cleaned up.
