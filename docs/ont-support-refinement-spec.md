# Exact early rejection experiment

The indexed full-input pass still processes roughly 100 batches per two minutes
on the large ONT input. Keep that producing checkout frozen while evaluating a
separate refinement. The old span acceptance rule requires every shared block
to have the same read-minus-path query offset: equal intervening lengths imply
equal offsets. Reject immediately when any two shared blocks disagree, instead
of collecting all remaining intersections before failing that same condition.
Preserve the original count, order and unplaced-sequence checks for survivors.

Compare with the preceding matcher on randomized geometry, regression fixtures
and frozen real candidate paths/read observations. Measure exact support work
separately; do not infer a complete reconstruction speedup from the experiment.
Only adopt a verified improvement in a separately declared producing attempt.

The real batch profile also identifies per-base placement comparisons in
`_witnesses`. Compare overlapping query-coordinate runs instead: within a run,
contig, strand and the genomic position at its first overlap determine every
placed coordinate. Ignore an interval where either side is unplaced, exactly
as the existing per-base predicate does. Keep exact DNA and anchor checks.
Validate this run comparison against the literal old predicate on random
placements, shifts and intervals, then compare complete witness lists and the
whole consumed-batch state on real paths and original reads.

Merge candidate genomic blocks per contig for the initial overlap screen, while
still testing each actual CIGAR block so skipped reference gaps remain skipped.
Bound sequence-match and interval memoization at 4,096 entries/offsets; cache
exhaustion must fall back to the exact calculation and retain every occurrence.
Hash benchmark results incrementally to avoid a second full serialized copy.
These refinements address [#486](https://github.com/openvax/isovar/issues/486).
