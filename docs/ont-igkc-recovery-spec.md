# Complete the two unresolved ONT views (#470)

Use the streaming engine from #483 in a separately identified offline attempt.
Verify the exact preserved BAM/index, acquisition receipt, source identity,
inventory, transcript references and both oriented reconstruction inputs before
running. Keep discovery/scientific settings unchanged and pin the new engine.
Preserve the original timeout and interrupted outcome with unknown counts.

Run both orientations serially with a declared one-hour reconstruction budget
per view, a 2 GiB scratch limit, a 16 GiB peak-RSS limit and an 8 GiB free
disk reserve. Write controls before computation, checkpoint each outcome, and
record time, peak memory, scratch use and the original incomplete-partner flags.
Earlier diagnostic attempts remain separate from this producing attempt.

Independently recount original records and validate candidate translation,
sequence witnesses and support identities before changing any coverage outcome.
Compare frozen discovery candidate sequences against the preceding engine and
report explicitly where the original failed attempt has no comparable result.

Publish a new selected-recovery ledger, revised pilot and event report retaining
all previous outcomes and producing requests. Preserve mixed engine provenance
for the unchanged 42 completed views and the two new reconstructions. Update
portable replay plans and validate a packaged offline replay. Do not convert
timeouts, missing partners or bounded discovery into negative biological claims.

Run lint, full tests and CI, then merge and release the separate recovery PR.

## Re-plan after the first streaming attempt

Preserve the operator-interrupted first attempt and its corrected diagnosis.
The initial missing-run-cache hypothesis was incorrect: coordinate objects
already cache their runs. CPU sampling instead identified repeated span work;
the current span matcher compares every read block with every candidate block
on the same contig/strand, including disjoint intervals. Index frozen blocks by
start and prefix maximum end, retaining nested/repeated placements and the exact
existing interval tests. Validate against the preceding brute-force matcher and
whole-result fixtures before declaring a new producing attempt. Keep the first
attempt's unknown support and all earlier failures in its provenance.

## Producing attempt after exact query and row sharing

Preserve the index-only forward one-hour timeout, its operator-interrupted
reverse comparison, and the separately recorded interruption of the preceding
refinement. After full validation, declare a new producer with the same 14,400
seconds of reconstruction time used by the original failed attempt. Keep the
2 GiB scratch, 16 GiB peak-RSS and 8 GiB disk-reserve limits, extending
resource monitoring through result export (#487). Partial exports remain
explicit failures. Keep scientific settings and original input pins unchanged.
Record monotonic reconstruction time and CPU time; host suspension is not
included in the monotonic timings on this macOS host.

## Retry after streaming audit publication

Preserve the subsequent forward free-space failure and the fact that reverse
never started. After #488 ships, merge that release into a new producer checkout
and freeze a separate environment from its verified dependency artifacts. Recheck
source, acquisition, reference, BAM/index and both orientation identities, plus
an independent native input count, before creating new controls. Keep the same
scientific settings and the 14,400-second, 2 GiB scratch, 16 GiB peak-RSS
and 8 GiB free-space budgets. Keep the machine awake for the serial runs.

Pin and retain the setup/coordinator scripts with the producing controls. Check
every completed candidate against original SAM evidence and independent ORF
translation. The preceding engine provides a baseline for all assembled path
IDs, DNA and blocks; it has no completed ONT full-support or ORF result. Report
that limitation explicitly rather than treating either prior failure as zero.

Measure output size and parsing memory before publication and installed replay:
streamed hashing/writing does not make the standard JSON decoder bounded in
memory. Preserve the 42 completed views and their producing identities, add the
two verified successors, and update coverage, reports and portable replay only
after both orientations pass. Retain the four existing archives and publish the
new evidence as a separately verifiable addition with an actual restore check.

## Resource measurement scope

The producing driver measures `resource.getrusage(...).ru_maxrss`, so the
16 GiB guard is a peak resident-set (RSS) limit. Report total process footprint
and compressed/swapped memory separately when available; the RSS limit does
not bound those quantities. Reconstruction, fingerprinting, serialization and
restoration also need separate timings and storage measurements.

## Validation and archive staging after incremental decoding

Keep the producing checkout, wheel, dependency environment and controls frozen.
Use the separately pinned incremental reader from #497 for independent checks;
record its identity alongside the producing engine. Validate every original SAM
record and native ORF before importing either successor. Run the checks serially
after both producers exit, under separately recorded validation resource limits.
The validation projection checks each complete path and every native ORF before
retaining only the fields used by official event-ORF export. Verify its complete
summary and original-record hashes against materialized fixtures before use.
The independent checker additionally limits sampled physical footprint to
24 GiB, because resident memory alone excludes compressed/swapped allocations.

Event reports should read one path at a time and retain only its geometry.
The recovery report should read only the observations and original SAM rows
used by its selected candidate, using separate passes so JSON field order does
not change the result. Validate all JSON, including skipped fields; preserve
every source file checksum and sequence/support check.

The publication and installed replay wheel includes #497 and #499. Record that wheel's
identity separately from the original producing wheel and verify complete result
equivalence on replay. Preserve the original producer request in the bundle.

Stage the two orientations separately for archival. Each addition must retain
its native result, normalized reconstruction, original input, producing controls,
previous outcomes and checks. Shared reports must identify their dependencies
on the preserved archives. Measure the actual payload before each allocation;
perform an actual archive restore and checksum verification before reclaiming
any session-created duplicate. Keep the four earlier archives unchanged.
Retain all unique evidence in verified archives and record any working-copy
cleanup explicitly. Complete installed full-input replay and checkpoint resume
before removing their disposable outputs.

The final replay writes a deterministic compressed checkpoint. Its separate CLI
resume must stream and validate the entire result checksum and comparison;
record unchanged checkpoint bytes before reclaiming that disposable output.
Preserve native and normalized evidence in one archive per orientation, plus a
small final report/provenance archive that resolves every combined report input
across the six data archives. Restore each archive through the public API before
reclaiming working duplicates. For reporting, restore only the normalized files
needed from those verified archives. Test exact-member restoration and cleanup
rejection on corrupt/changed fixtures before using them on producing evidence.

Freeze the final installed replay environment separately from the producer,
including its updated repository release dependencies. Compare every installed
package module with the wheel and source and explicitly confirm that changes
from the producer are confined to audit serialization and verification. Preserve
the producer environment and historical dependency/reference identities.
