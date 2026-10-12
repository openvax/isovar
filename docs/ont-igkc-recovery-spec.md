# Complete the two unresolved ONT views (#470)

Use the streaming engine from #483 in a separately identified offline attempt.
Verify the exact preserved BAM/index, acquisition receipt, source identity,
inventory, transcript references and both oriented reconstruction inputs before
running. Keep discovery/scientific settings unchanged and pin the new engine.
Preserve the original timeout and interrupted outcome with unknown counts.

Run both orientations serially with a declared one-hour reconstruction budget
per view, a 2 GiB scratch limit, a 16 GiB process-memory limit and an 8 GiB free
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
2 GiB scratch, 16 GiB process-memory and 8 GiB disk-reserve limits, extending
resource monitoring through result export (#487). Partial exports remain
explicit failures. Keep scientific settings and original input pins unchanged.
Record monotonic reconstruction time and CPU time; host suspension is not
included in the monotonic timings on this macOS host.
