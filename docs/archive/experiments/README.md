# Historical experiments

These reports preserve the inputs, measurements and code decisions behind
earlier Isovar changes. They record what a specific run established and its
limits; current options and interpretation belong in the main documentation.

| Experiment | Measured conclusion and code change |
|---|---|
| [Dense SV evidence-group budget (#464)](sv-evidence-groups-464.md) | A group cap bounded discovery resources on selected COLO829/HCC1395 RNA, while complete original-read support remained separate. Led to the optional `max_evidence_groups` limit and explicit omission counts in PR #567. |
| [SV span and witness queries (#485)](sv-span-and-witness-queries-485.md) | An interval index eliminated disjoint candidate-block scans in a component test; later exact-query sharing preserved full-result equality. Led to `_PathSpans` indexing and bounded witness-query caches. |
| [High-depth RNA support matching (#227)](osteosarc-rna-performance/README.md) | T3 deduplicated MT-ND5 fell from 1,386.0 to 69.2 seconds with identical ordered outputs in completed comparisons. Led to indexed support matching in PR #230; host timings are not a portable guarantee. |
| [Read-collection memory (#229)](osteosarc-rna-performance/COLLECTION_MEMORY.md) | Deep collection peaks fell about 70% with unchanged evidence and higher CPU. Led to per-call coordinate sharing in PR #232, later refined in PR #239. |
| [Coordinate sharing (#234)](osteosarc-rna-performance/COORDINATE_SHARING.md) | The old LRU bound was not reached in the scanned real loci; a dict kept memory within 3% and reduced collection CPU about 9% across nine pairs. Led to PR #239; no whole-pipeline speedup was claimed. |
