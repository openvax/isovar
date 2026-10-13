# Dense reconstruction performance

Follow-up to #464 / #475, implemented by #481. The target is a 10x runtime
improvement. Report measured reconstruction and output costs separately; neither
compression ratios nor a timed-out baseline establish an end-to-end speedup.

## Method

`examples/benchmark_dense_sv.py` runs one preserved audit-bundle input in a fresh
process, verifies its BAM/index/input checksums before timing, disables network
access, and records reconstruction wall/CPU time, peak RSS, sampled scratch,
discovery selection and per-candidate sequence/support fingerprints. Its optional
profiler excludes output serialization. Original input-acquisition limitations
remain explicit. Compact gzip result serialization is timed separately and is
not the CLI's uncompressed JSON-writing cost.

The baseline is released 1.49.0 (`546b95a`, the engine merged as `21b5020`). Both
versions use the same repository-local Python 3.11 environment, inputs and
existing discovery settings. These are single runs on a shared macOS arm64
machine, not confidence intervals. Required local tests overlapped the final
current-engine runs; an independent diagnostic overlapped the 100k baseline.
CPU time and raw timing records accompany wall time.

The inputs and pins were previously published in [PR #475](https://github.com/openvax/isovar/pull/475).
The three pinned views come from bundle
`744b38a81f2c8775896bc21e5975bea148e006ec10051ac915afe2c8528e8630`:

- Short reads: source `6e9147690525529f9591602afd1191f9862f4df16117a865c6726654456c9368`.
- ONT: source `58f81602ec2f87b86a032cf869f1abefc9f4235dbb1ce082e0b6bbc085b155de`.
- Smaller forward view: `adj-03b9fc802bfa90900b74a750` for both sources.
- Capped short-read IGKC forward view: `adj-041dbfbb2239d183a5982838`.

The small views fit inside the 50,000-record discovery limit; IGKC exercises that
limit. The same original voters, budget semantics and full-input support must
survive optimization. Grouping is computational reuse, not molecular
or biological deduplication.

## Remaining full-input support work

The large recovery remains unresolved in #470. A separate diagnostic did not
produce a completed result, so it establishes no end-to-end speedup. Its private
input identifiers and diagnostic measurements remain outside this PR.

The next support change (#482) must release transient per-base observations after
recording their contributions. The current pass retains qualifying templates
and their original/mirrored geometries before reporting complete-input support;
this still scales with total qualifying aligned bases.

## Correctness checks

- Full lint and test-suite gates, including exact CIGAR geometry, original mate
  labels, missing qualities, alternative placements, late rare events, input
  immutability, SAM/BAM payload retention and deterministic overlap caps.
- 2,500 randomized CIGAR/strand/quality/clip cases matched released geometry.
- 40 randomized complete reconstructions matched released extension, including
  alternate placements and reverse paths.
- All four complete result objects match 1.49.0, including observations, original
  records, qualities/metadata evidence, paths, discovery and support summaries.
- Final local gates: `./lint.sh`, `./test.sh` (5,899 passed), and `pip check`.

## Measured reconstruction results

The 10x overall target was **not reached**. Gains vary with input and discovery
budget; these measurements do not establish a general throughput guarantee.

| View | Record budget | 1.49.0 wall / CPU (s) | Current wall / CPU (s) | Wall speedup | Complete result |
|---|---:|---:|---:|---:|---|
| Short reads | 50,000 | 25.08 / 25.04 | 21.98 / 22.02 | 1.14x | Identical |
| ONT | 50,000 | 28.67 / 28.63 | 18.57 / 18.59 | 1.54x | Identical |
| IGKC short reads | 50,000 | 192.29 / 186.27 | 59.68 / 59.77 | 3.22x | Identical |
| IGKC short reads | 100,000 | 360.50 / 348.87 | 83.28 / 83.41 | 4.33x | Identical |

| View / budget | 1.49.0 RSS (MiB) | Current RSS (MiB) | 1.49.0 scratch (MiB) | Current scratch (MiB) |
|---|---:|---:|---:|---:|
| Short reads / 50k | 605.9 | 634.2 | 46.3 | 16.5 |
| ONT / 50k | 647.5 | 696.5 | 53.7 | 21.7 |
| IGKC short reads / 50k | 2467.7 | 2454.6 | 785.0 | 242.3 |
| IGKC short reads / 100k | 2398.1 | 2706.7 | 788.8 | 246.2 |

Scratch use falls substantially. Peak RSS is roughly unchanged to 13% higher in
these single runs; this PR does **not** establish a general memory reduction.
The streaming support follow-up is [#482](https://github.com/openvax/isovar/issues/482).
Discovery remains bounded and complete-input support is scoped to the preserved
input, whose upstream acquisition can itself be incomplete.

[Machine-readable measurements](dense-group-execution-benchmark.json) include
input and engine hashes, parameters, stage costs, output-serialization cost and
whole-result/candidate-support comparison digests. Full result artifacts remain
in the local benchmark directory; no original read payloads are added here.

To repeat a case with the same driver and its report's run ID:

```sh
python examples/benchmark_dense_sv.py --checkout /path/to/checkout \
  --bundle /path/to/preserved-bundle --run 'SOURCE/ADJACENCY/forward' \
  --max-records 50000 --seconds 900 --max-scratch-bytes 8589934592 \
  --output result.json
```

Repeat with the baseline checkout, using identical input and parameters. The
driver writes both measurements and a compressed complete result. Compare full
result objects as well as the candidate/support fingerprints; a capped discovery
search must not be described as exhaustive.
