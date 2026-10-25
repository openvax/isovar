# Historical experiments

These reports preserve the inputs, measurements and code decisions behind
earlier Isovar changes. They record what a specific run established and its
limits; current options and interpretation belong in the main documentation.

| Experiment | Measured conclusion and code change |
|---|---|
| [Dense SV evidence-group budget (#464)](sv-evidence-groups-464.md) | A group cap bounded discovery resources on selected COLO829/HCC1395 RNA, while complete original-read support remained separate. Led to the optional `max_evidence_groups` limit and explicit omission counts in PR #567. |
| [SV span and witness queries (#485)](sv-span-and-witness-queries-485.md) | An interval index eliminated disjoint candidate-block scans in a component test; later exact-query sharing preserved full-result equality. Led to `_PathSpans` indexing and bounded witness-query caches. |
