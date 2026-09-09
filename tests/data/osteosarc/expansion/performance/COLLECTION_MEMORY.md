# Read-collection memory (#229): implementation spec

## Scope

Reduce high-depth RNA collection peak memory without changing biological
semantics, evidence units, filtering, deletion normalization, mate merging,
read ordering, or original read identities/multiplicity. Preserve public
list-returning methods and existing collector subclass hooks. This follows
the reconstruction optimization released in Isovar 1.8.2 / PR #230.

## Plan and acceptance

1. Profile allocations separately from uninstrumented timing/RSS on the
   checksum-verified full-region original RNA sources. Attribute retained
   allocations before choosing a storage/lifetime change.
2. Minimize retained per-base collection data. Keep the existing conversion
   and merging semantics; document any new internal representation or
   extension-hook compatibility handling.
3. Differential tests against frozen 1.8.2 collection cover exact ordered
   locus/allele reads and ref/alt/other evidence, conflicting/overlapping mates,
   supplementary alignments, shared names, quality/flag/soft-clip policies,
   indels and empty results. Add deterministic lifetime/storage regressions.
4. Compare fresh-process full-depth collection and complete RNA/protein
   pipeline timing, peak RSS and exact output fingerprints against 1.8.2.
   Include original and primary-only modes, long-read cases and a bulk control.
   Setup and profiling overhead must be labeled, not hidden as improvements.
5. Bump to 1.8.3, run lint/full tests and CI, and open a PR for review.
   Do not merge or publish as part of this PR-start request.

The full source cache is acquisition input, not a new checked-in BAM corpus.
Preserve the existing historical benchmark/audit artifacts unchanged and
record new measurements separately.

## Allocation attribution and selected design

The initial full-depth T1 tagged allocation trace on unchanged 1.8.2 retained
1,298,324,063 traced bytes before mate merging. Coordinate integers allocated
by `get_reference_positions(full_length=True)` accounted for 977,729,120 bytes
(30,554,035 allocations), with position-list slices another 246,861,808 bytes.
This instrumented run used 6.77 GB peak RSS including profiler/snapshot overhead;
it is not an uninstrumented memory or runtime baseline.

Share equal built-in coordinate integers in a per-`get_locus_reads` bounded
LRU cache (65,536 entries). Each read retains its own existing list object;
custom coordinate objects and sequence types returned by subclass hooks are
left untouched. Cache eviction affects sharing only, never values or evidence.
The cache is local to one call, not persistent collector state. No streaming
API, compressed coordinate type, hook bypass, or mate-merging change is needed.
