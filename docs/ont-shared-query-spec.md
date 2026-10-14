# Share exact interval searches within each support batch

The forward ONT discovery paths generate 899 interval queries but only 84
unique DNA strings. Search each unique interval against the batch's distinct
observed sequences once, then reuse the exact offsets across candidate paths.
Discard nonmatching observations before the unchanged anchor and placement
checks. Preserve each junction's original observation membership, all repeated
matches, original identities and witness ordering.

Keep the index batch-local and cap retained queries, sequences and offsets at
4,096 entries. Oversized/repetitive results or a full cache fall back to the
existing exact search, with no occurrence limit. Do not alter either running
producer checkout. Compare complete witness lists, full consumed-batch state,
whole-result fixtures and the public short-read benchmark before considering a
separately identified producing attempt.

Check placed coordinates within the exact DNA interval before searching for
flank anchors. A mismatch within that interval necessarily fails the existing
whole-link placement rule. Reuse this result when the eventual link is the same
interval; extended links still receive the full sequence/placement check.

A later input batch also produced 2,326 witness rows but only 460 distinct row
values. Intern complete equal witness rows within a batch, capped at 4,096
unique rows, including original observation identity and every junction link in
the key. Retain shared rows when the final junction filter removes no links.
This preserves every output value/count; it only shares equal Python objects.
