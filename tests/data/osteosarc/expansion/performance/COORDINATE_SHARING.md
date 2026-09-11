# Per-call coordinate sharing (#234)

## Change

Isovar 1.8.3 (#229, #232) shared equal coordinate integers across the reads of
one `ReadCollector.get_locus_reads` call through a `functools.lru_cache`
bounded at 65,536 entries. Reads arrive in coordinate order and each read's
positions increase, so once a locus has more distinct coordinates than the
bound, every lookup evicts an entry the next read needs. Sharing then drops to
zero while every lookup still pays LRU bookkeeping.

Isovar 1.8.4 shares through a plain dict local to the call. Every cached
integer is already held by a read's coordinate list, so the dict adds only
table slots, and its size is bounded by the distinct aligned coordinates of
the reads at that locus rather than by a fixed count. Coordinate values, order
and list objects, hook calls, mate merging and evidence are unchanged. Only
the mechanism that finds the shared integer differs.

## Corpus bound scan

Before measuring, every nuclear inventory variant was collected from every
source in the checksum-verified region cache with unmodified 1.8.3, with its
LRU cache instrumented. None of the 4,563 non-empty nuclear source/locus pairs
exceeds the bound. The largest working set, EXOC4 chr7:133274996 in the T3
tagged source, has 46,702 distinct coordinates, and the LRU served 98.6% of
lookups there from cache. Mitochondrial loci cannot reach the bound, because
the whole mitochondrial genome has 16,569 positions.

The recorded corpus therefore does not exhibit the zero-sharing case. The
regression test `test_reads_longer_than_the_former_lru_bound_still_share_coordinates`
reproduces it deterministically. The full-depth runs below measure the effect
of the new mechanism on real data, including the near-bound EXOC4 locus. The
scan is a diagnostic for choosing loci, not packaged evidence.

## Results

All ten fresh-process collection comparisons against 1.8.3 preserve every
field in the ordered locus-read and allele-read fingerprints, including
coordinate and quality container types, sequence, intervals, names and source
counts, and all ref/alt/other counts match. RSS is the process high-water at
return from `get_locus_reads`, before fingerprinting and conversion. GB means
decimal GB.

| Full source / locus | Returned locus objects | 1.8.3 peak GB | 1.8.4 peak GB | 1.8.3 CPU s | 1.8.4 CPU s | CPU change |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| T1 tagged / MT-ND5 | 34,026 | 0.490 | 0.492 | 4.874 | 4.601 | −5.6% |
| T1 tagged, primary-only / MT-ND5 | 33,800 | 0.487 | 0.485 | 4.755 | 4.453 | −6.3% |
| T3 deduplicated / MT-ND5 | 105,522 | 1.483 | 1.491 | 20.292 | 17.244 | −15.0% |
| T3 deduplicated, primary-only / MT-ND5 | 105,183 | 1.489 | 1.488 | 18.539 | 18.596 | +0.3% |
| T3 tagged / MT-ND5 | 199,767 | 2.697 | 2.695 | 31.445 | 28.094 | −10.7% |
| T1 tagged / DYNC1H1 chr14:101980529 | 5,007 | 0.228 | 0.227 | 0.908 | 0.854 | −5.9% |
| T2 tagged / H1-2 deletion | 21,439 | 0.323 | 0.322 | 2.077 | 1.923 | −7.4% |
| Bulk 2024 / DYNC1H1 chr14:102030200 | 452 | 0.146 | 0.147 | 0.025 | 0.024 | −3.6% |
| T3 tagged / EXOC4 chr7:133274996 | 2,370 | 0.194 | 0.189 | 0.453 | 0.388 | −14.4% |

Peak memory is unchanged: every pair is within 3%, as expected where the LRU
already shared nearly every coordinate. Collection CPU falls on eight of nine
pairs, by 4–15%, and is flat on T3 deduplicated primary-only; across all nine
pairs it falls by about 9%. These are single observations per side, and the
whole-pipeline collection stage below shows no clear difference, so this is not
a speedup claim. Isovar 1.8.3 raised collection CPU by 36–70% over 1.8.2, and
1.8.4 removes only the LRU bookkeeping, not the per-read list rebuild. Wall
times on this shared macOS ARM64 / Python 3.12.6 host vary with load: the T3
deduplicated primary-only 1.8.4 run took 33% more wall time than its baseline
at equal CPU.

The allocation-instrumented T1 pair traces the same 34,071 pre-merge reads and
30,619,733 query bases on both sides. Retained traced memory before mate
merging is 322,395,855 bytes on 1.8.3 and 321,621,791 bytes on 1.8.4, 0.24%
lower, so the dict adds no retained memory relative to the bounded LRU. As in
#229, instrumented runs are allocation attribution, not timing or RSS
measurements.

Two whole-pipeline T3 tagged MT-ND5 pairs, the second with the run order
reversed, preserve exact ordered RNA candidates, coverage and original-read
fingerprints, and all 4,351 independently checked ranked protein windows. Peak
RSS is the process high-water through independent validation, as in #229.

| Pair | 1.8.3 peak GB | 1.8.4 peak GB | 1.8.3 production CPU s | 1.8.4 production CPU s | 1.8.3 check CPU s | 1.8.4 check CPU s |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| First, 1.8.3 ran first | 3.523 | 3.868 | 86.3 | 86.5 | 13.7 | 11.9 |
| Repeat, 1.8.4 ran first | 4.773 | 3.633 | 87.1 | 84.2 | 18.1 | 14.6 |

Whole-pipeline peak RSS on this host is dominated by host memory state, not by
this change. The run order decides which side peaks higher, and the identical
1.8.3 code and inputs peaked at 3.523 GB and 4.773 GB in this session and at
6.173 GB in the #229 session. The first 1.8.4 pipeline run was also contended:
its read-collection stage used 45.6 s of wall time for 38.1 s of CPU. In the
uncontended repeat, that stage matches closely, peaking at 3.183 GB with
32.4 s of CPU on 1.8.3 and 3.177 GB with 32.0 s on 1.8.4. Collection returns
identical reads, and the standalone T3 tagged benchmark's whole-process peak
after conversion is 3.206 GB on 1.8.3 and 3.204 GB on 1.8.4. No whole-pipeline
memory or CPU change is claimed in either direction.

## Evidence

[`coordinate-sharing-results/manifest.json`](coordinate-sharing-results/manifest.json)
pins [`collection.json`](coordinate-sharing-results/collection.json),
[`pipeline.json`](coordinate-sharing-results/pipeline.json) and the digest of
the report generator that packaged them. CI verifies the package, the exact
before/after fingerprints and the frozen observations above without needing
full-depth inputs or imposing machine-dependent timing gates.

Baseline runs import a clean detached `v1.8.3` worktree and final runs import
this branch. Both sides of every pair run the same `collection_benchmark.py`
from this branch, so each pair shares one benchmark digest. That script now
records sequence container types in its read fingerprints (#237), so these
runs cannot be paired with the #229 records, and the 1.8.3 baselines were
measured afresh on the same host in the same session. Each pair alternates
which side runs first. Every run hashes its full input BAM before collection,
so both sides start with the same file cache state.

## Reproduction

Obtain the checksum-verified original region cache as described in the parent
expansion README, and create a clean detached `v1.8.3` worktree. Run from a
neutral working directory with `PYTHONPATH` pointing at the intended checkout,
exactly as in [`COLLECTION_MEMORY.md`](COLLECTION_MEMORY.md), always invoking
this branch's `collection_benchmark.py` by path:

```sh
PYTHONPATH="$baseline_worktree" python \
  "$repo/tests/data/osteosarc/expansion/collection_benchmark.py" \
  --destination "$sources" --source-id 58c4d68e5f7e55b5 \
  --variant-id EXOC4-chr7-133274996 --output "$runs/baseline-t3exoc4"
PYTHONPATH="$repo" python \
  "$repo/tests/data/osteosarc/expansion/collection_benchmark.py" \
  --destination "$sources" --source-id 58c4d68e5f7e55b5 \
  --variant-id EXOC4-chr7-133274996 --output "$runs/final-t3exoc4"
```

| Label | Source | Variant | Mode |
| --- | --- | --- | --- |
| t1 | 53f498a544883d51 | MT_ND5-chrM-12994 | defaults |
| t1primary | 53f498a544883d51 | MT_ND5-chrM-12994 | primary_only |
| t3dedup | 62a9d3e67e1b564e | MT_ND5-chrM-12994 | defaults |
| t3dedupprimary | 62a9d3e67e1b564e | MT_ND5-chrM-12994 | primary_only |
| t3tagged | 58c4d68e5f7e55b5 | MT_ND5-chrM-12994 | defaults |
| ontdync | 53f498a544883d51 | DYNC1H1-chr14-101980529 | defaults |
| ontdel | 1120a096937e29e1 | H1_2-chr6-26055824 | defaults |
| bulk | 1b66c15da594a3ef | DYNC1H1-chr14-102030200 | defaults |
| t3exoc4 | 58c4d68e5f7e55b5 | EXOC4-chr7-133274996 | defaults |
| profile_t1 | 53f498a544883d51 | MT_ND5-chrM-12994 | defaults, `--allocations` |

The whole-pipeline pair uses `python -m tests.data.osteosarc.expansion.benchmark`
on the T3 tagged MT-ND5 source with `--timeout 600`, importing each checkout's
own unchanged copy of that module. Package with `python -m
tests.data.osteosarc.expansion.collection_report`, passing each run as
`--collection baseline_<label>=<dir>` or `--collection final_<label>=<dir>`,
the pipeline runs as `--pipeline baseline_t3tagged=<dir>` and
`--pipeline final_t3tagged=<dir>`, and each pair as a `--compare-collection` or
`--compare-pipeline` argument. `--output` must be a new directory.
