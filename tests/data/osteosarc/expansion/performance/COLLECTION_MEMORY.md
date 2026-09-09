# Read-collection memory (#229)

## Uninstrumented collection results

Eight fresh-process comparisons against unchanged 1.8.2 preserve every field
in the ordered locus-read and allele-read fingerprints, including coordinate
lists, qualities, sequence, query/reference intervals, names and source counts.
All ref/alt/other counts match. These are full-region original inputs, not the
small checked-in fixtures. RSS below is the process high-water at return from
`get_locus_reads`, before fingerprinting/conversion. GB means decimal GB.

| Full source / locus | Returned locus objects | 1.8.2 peak GB | New peak GB | 1.8.2 CPU s | New CPU s |
| --- | ---: | ---: | ---: | ---: | ---: |
| T1 tagged / MT-ND5 | 34,026 | 1.466 | 0.480 | 4.107 | 5.606 |
| T1 tagged, primary-only / MT-ND5 | 33,800 | 1.483 | 0.506 | 2.969 | 4.440 |
| T3 deduplicated / MT-ND5 | 105,522 | 5.173 | 1.513 | 11.680 | 17.295 |
| T3 deduplicated, primary-only / MT-ND5 | 105,183 | 5.163 | 1.514 | 11.366 | 16.988 |
| T3 tagged / MT-ND5 | 199,767 | 9.448 | 2.692 | 23.643 | 33.875 |
| T1 tagged / DYNC1H1 chr14:101980529 | 5,007 | 0.440 | 0.222 | 0.604 | 1.024 |
| T2 tagged / H1-2 deletion | 21,439 | 0.809 | 0.318 | 1.599 | 2.496 |
| Bulk 2024 / DYNC1H1 chr14:102030200 | 452 | 0.143 | 0.140 | 0.022 | 0.026 |

This is a memory/CPU tradeoff, not a speedup. The deepest collection peaks
fall by about 70%, while collection CPU rises roughly 43–49% on these deep
inputs. Absolute bulk differences are tiny. Wall times on this shared macOS
ARM64/Python 3.12.6 host vary with load; the full records retain both wall and
CPU time, plus post-collection work and whole-process peak RSS. Allocation-
instrumented runs must not be substituted for these ordinary measurements.

Read objects and raw alignment/template counts are different evidence units,
not cell or UMI counts. No coordinate cache entry is a vote or observation.

## Full RNA/protein pipeline results

All seven fresh-process comparisons pass exact ordered RNA candidate/coverage/
original-read fingerprints and complete ranked protein/check/support-name
fingerprints. Each side independently checks all **9,948 ranked windows** in
total. The ABCF2 control correctly records zero candidates and checked proteins;
its underlying read counts are not zero and remain identical.

| Full source / locus | 1.8.2 peak GB | New peak GB | 1.8.2 production CPU s | New production CPU s | New independent-check CPU s | Ranked windows |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| T1 tagged, primary-only / MT-ND5 | 1.606 | 1.124 | 9.295 | 10.600 | 1.549 | 2,006 |
| T3 deduplicated / MT-ND5 | 5.496 | 3.061 | 30.485 | 34.975 | 5.572 | 2,854 |
| T3 tagged / MT-ND5 | 9.954 | 6.173 | 75.786 | 76.305 | 13.928 | 4,351 |
| T1 tagged / DYNC1H1 chr14:101980529 | 0.508 | 0.292 | 0.750 | 1.053 | 0.019 | 241 |
| T2 tagged / H1-2 deletion | 0.897 | 0.423 | 1.963 | 2.613 | 0.039 | 384 |
| Bulk 2024 / DYNC1H1 chr14:102030200 | 0.198 | 0.196 | 0.040 | 0.045 | 0.002 | 112 |
| Bulk / ABCF2 (no RNA candidates) | 0.194 | 0.195 | 0.014 | 0.016 | <0.001 | 0 |

Peak RSS here is the process high-water through independent validation,
including earlier setup/collection/reconstruction and diagnostic RNA hashing;
it excludes the later compressed protein-report serialization. Timings exclude
reference setup, checksum verification, RNA hashing and report serialization.
Production and independent checking are measured separately in the JSON.
The deep whole-pipeline peaks fall by **30–44%**, not 70%: after the collection
fix, reconstruction/checking becomes the new high-water. Production CPU rises
about 14–15% on T1-primary/T3-deduplicated; T3-tagged is approximately unchanged
in this observation. Nuclear long-read CPU also rises; this PR makes no
whole-pipeline speedup claim.

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

The recorded full-depth T1 tagged allocation trace on unchanged 1.8.2 retained
1,298,324,919 traced bytes before mate merging. Coordinate integers allocated
by `get_reference_positions(full_length=True)` accounted for 977,729,120 bytes
(30,554,035 allocations), with position-list slices another 246,862,032 bytes.
The optimized trace retains 322,395,887 bytes, including just 436,160 bytes of
those coordinate allocations. Both contain 34,071 pre-merge reads and
30,619,733 query bases, with identical post-merge fingerprints and counts.

Tracer bookkeeping alone is 1,524,085,808 bytes before and 18,754,160 after.
Process high-water including tracing/snapshot analysis is 8.51 GB before and
0.59 GB after; these are **not** ordinary memory/runtime measurements. Tracing
is stopped before snapshot aggregation, and error paths restore tracer state.

Share equal built-in coordinate integers in a per-`get_locus_reads` bounded
LRU cache (65,536 entries). Each read retains its own existing list object;
custom coordinate objects and sequence types returned by subclass hooks are
left untouched. Cache eviction affects sharing only, never values or evidence.
The cache is local to one call, not persistent collector state. No streaming
API, compressed coordinate type, hook bypass, or mate-merging change is needed.

## Validation follow-up (#233)

Full-suite scheduling exposed an existing fixture-identity collision: the
24-transcript stress subset and 131-transcript cohort declared the same
assembly/release dataset label. Varcode's name-keyed contig cache
(openvax/varcode#402) could then reject valid chr6/chr7 queries depending on
which dataset ran first. The offline reference loader now appends the pinned
manifest digest to its runtime dataset name. Original assets/manifests remain
unchanged; both loading orders and genuine gene lookups are tested. This
repairs Isovar's ambiguous fixture identity, not Varcode's general cache.

## Evidence and interrupted attempts

[`collection-results/manifest.json`](collection-results/manifest.json) pins
[`collection.json`](collection-results/collection.json) and
[`pipeline.json`](collection-results/pipeline.json). They retain all completed
measurements, identities and raw-artifact checksums; the full protein dumps
and original BAMs remain outside git. The completed comparisons comprise
eight ordinary collection pairs, one explicitly instrumented pair, and seven
full-pipeline pairs. CI verifies the packaged evidence and checksum manifest
without needing full-depth inputs or imposing machine-dependent timing gates.

The report also preserves two unsuccessful exploratory runs, excluded from
successful comparisons: an interrupted original allocation profile, and a
pipeline run whose compressed diagnostic dump failed with `ENOSPC`. The latter
also could not write its final result record; its marked error-only record was
recovered from the captured exception log after disk space became available.
It has no inferred counts or completion fingerprint. Its paired initial
baseline remains as an unpaired observation. No source or partial output was
deleted; fresh directories were used for all successful retries. A full test
attempt during the disk-full condition was likewise discarded and rerun.

## Reproduction

Use the acquisition instructions in the parent expansion README to obtain
the checksum-verified original region cache. No credentials or full-depth
BAMs are committed. The recorded identities contain every source/locus/mode,
input BAM/index/receipt/inventory digest, interpreter/platform and imported
source-file digest. These are single fresh-process observations on one host,
not a cross-platform performance guarantee or a statistical benchmark.

Create a clean detached `v1.8.2` worktree. Run from a neutral working directory
(not either checkout), with `PYTHONPATH` pointing to the intended checkout so
an editable installation cannot accidentally supply the other implementation.
For example, with absolute paths in `collection_repo`, `collection_baseline`,
`collection_sources`, and `collection_output`:

```sh
PYTHONPATH="$collection_baseline" python \
  "$collection_repo/tests/data/osteosarc/expansion/collection_benchmark.py" \
  --destination "$collection_sources" --source-id 58c4d68e5f7e55b5 \
  --output "$collection_output/baseline-t3tagged"
PYTHONPATH="$collection_repo" python \
  "$collection_repo/tests/data/osteosarc/expansion/collection_benchmark.py" \
  --destination "$collection_sources" --source-id 58c4d68e5f7e55b5 \
  --output "$collection_output/final-t3tagged"
```

The default locus above is MT-ND5. Use `--variant-id` for nuclear controls and
`--mode primary_only` for the recorded alternative filtering mode. For each
full-pipeline pair, use `python -m
tests.data.osteosarc.expansion.benchmark` with the same source/locus/mode and
`--timeout 600`; that command checks every ranked window, not just the public
top result. Run allocation profiling separately with `collection_benchmark.py
--allocations` on T1 (`53f498a544883d51`). Do not run other heavy jobs alongside
benchmarks; retain interrupted outputs and use new directories for retries.

Package outputs with `python -m
tests.data.osteosarc.expansion.collection_report`: pass each collection run as
`--collection label=/absolute/run/directory`, each pipeline run as `--pipeline
label=/absolute/run/directory`, and each completed pair as
`--compare-collection baseline_label=final_label` or `--compare-pipeline
baseline_label=final_label`. `--output` must be a new directory. Comparisons
fail on changed evidence, inputs or incomplete runs; an error record must
never be paired as a successful empty result. Repackaging the same run
directories must reproduce the checked-in JSON and manifest byte for byte.
