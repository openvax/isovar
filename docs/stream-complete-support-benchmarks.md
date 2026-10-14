# Streaming complete-input support

PR #483 fixes #482 by consuming observations against frozen candidates in
batches of at most 250,000 observed bases (both orientations count). One
indivisible template can exceed this working-set bound. The existing geometry
cache remains bounded at two million coordinate entries; coordinate interning
also has a fixed bound. Neither bound changes discovery, evidence or support.
Original SAM payloads, observation rows, witnesses and consulted metadata remain
in the output and necessarily grow with supporting input.

## Matched public fixtures

Fresh processes compare released 1.52.0 with 1.53.0 on the same previously
published portable inputs used in PR #481. Each run has a 600-second runtime
budget, a 2 GiB sampled scratch limit and an 8 GiB sampled peak-RSS limit.
Times include reconstruction, excluding input verification and JSON serialization.
These are single shared-machine runs: parts overlapped local tests and synthetic
measurements, so timings are indicative rather than isolated speed estimates.

| Input / discovery cap | 1.52.0 time | Streaming time | 1.52.0 peak RSS | Streaming peak RSS |
|---|---:|---:|---:|---:|
| Smaller short-read view / 50k | 29.83 s | 26.27 s | 0.604 GiB | 0.572 GiB |
| Smaller ONT view / 50k | 19.83 s | 21.45 s | 0.681 GiB | 0.641 GiB |
| IGKC short reads / 50k | 67.90 s | 77.87 s | 2.402 GiB | 1.476 GiB |
| IGKC short reads / 100k | 87.38 s | 108.30 s | 2.644 GiB | 1.722 GiB |

All four **complete result objects match**, including every candidate, original
record, observation, witness, support count, label, parent and discovery field.
List order is preserved in the comparison. The two IGKC runs retain 189,190
observation rows while consuming 123 batches, none exceeding 250,000 observed
bases. Their peak RSS falls 35–39%, with a 15–24% runtime cost in these runs.
This is a memory improvement, not an end-to-end 10× speedup.

## Synthetic depth scaling and retained output

`examples/benchmark_stream_support.py` generates deterministic error-free reads
with varying flanks around one splice, freezes discovery at one record, and
counts support over every original. Each invocation has a 600-second budget.
The three baseline/current pairs have identical complete-result SHA-256 values.
The fixture deliberately includes different coverage lengths: every direct
junction read need not be an exact full-interval witness.

| Original reads | Baseline Python peak | Streaming Python peak | Retained output, baseline / streaming | Peak minus retained, baseline / streaming |
|---|---:|---:|---:|---:|
| 100 | 10.09 MB | 6.56 MB | 1.01 / 1.04 MB | 9.08 / 5.52 MB |
| 1,000 | 125.65 MB | 17.55 MB | 10.55 / 10.58 MB | 115.10 / 6.97 MB |
| 5,000 | 420.81 MB | 56.12 MB | 51.70 / 51.74 MB | 369.11 / 4.38 MB |

At 5,000 reads the Python allocation peak is 7.5× lower, with the same retained
output. Tracemalloc measures Python allocations during reconstruction, excludes
native libraries and serialization, and is stopped after garbage collection.
Peak minus retained is an estimate of excess working memory, **not** the peak
size of geometry alone: geometry and output can peak at different times. RSS
and serialization sizes are recorded separately. Weak-reference regression
tests also verify that streamed observation objects do not survive reporting.
These synthetic values are not measurements of the unresolved large ONT input.

The [machine-readable evidence](stream-complete-support-benchmark.json) contains
input pins, source/driver hashes, exact settings, stage/CPU/wall times, scratch,
RSS, batch occupancy, complete-result digests and the synthetic allocation data.

```sh
python examples/benchmark_dense_sv.py --checkout CHECKOUT --bundle BUNDLE \
  --run SOURCE/GEOMETRY/forward --max-records 50000 --seconds 600 \
  --max-scratch-bytes 2147483648 --max-rss-bytes 8589934592 --output measured.json
python examples/benchmark_stream_support.py --checkout CHECKOUT \
  --reads 5000 --seconds 600 --output synthetic.json
```

No original input payloads or private recovery diagnostic metadata are included
in this report. The separately pinned recovery and independent evidence checks
are documented in the [completed recovery](sid-sv-dense-recovery.md), #470 and
PR #484. These benchmark measurements retain their original inputs and versions.
