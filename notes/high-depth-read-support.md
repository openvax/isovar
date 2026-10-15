# High-depth read support (#463)

The original MT-ND5 workload spends time repeatedly sorting supporting reads.
A baseline cProfile run took 243 seconds, with 78 seconds inside sorting
(including its key functions). That profile is separate from the timing trials.

Two existing loops now avoid that repeated work:

- Shared-support assignment sorts each read group once, before visiting candidates.
- Exact-candidate merging first groups by sequence and transcript identity. It
  computes detailed read ordering only to break ties within those groups.

The same read ordering and alignment compatibility checks still apply. No reads
are capped or discarded, and there is no new API, cache, framework or dependency.

## Full-input comparison

Two fresh processes per implementation, with order reversed for the second pair,
using Python 3.14.6 on macOS and baseline commit
`c85823a11372dd0109435e95941ac89c55c29d09` (1.54.2):

| Measurement | Baseline | Changed |
| --- | ---: | ---: |
| Pipeline wall seconds | 105.2, 112.6 | 78.3, 69.8 |
| Pipeline CPU seconds | 103.7, 108.9 | 76.3, 69.2 |
| Peak process RSS, decimal GB | 1.57, 1.16 | 0.96, 1.01 |

Mean wall time fell from 108.9 to 74.1 seconds (32%). RSS varies with host memory
pressure; these are observed ranges, not a promised percentage. The benchmark
retains candidate and translation objects in both versions for later comparison.
Reference indexing and evidence fingerprinting are outside the measured pipeline.

All 1,851 candidates and 1,239 translations matched, including their order,
original read identities and metadata, chosen alignments, transcript constraints,
coverage, protein support and filters. Both versions retained 18,089 reads
(4,533 alternate) and 10,993 fragments (2,725 alternate). Small NME1, GTF3C5
deletion and AKT2 cases also matched completely.

The input is the original regional member from the published `openvax-v1` bundle,
not a new sample of its reads. Input, code and output hashes, dependency versions
and individual measurements are in
[the benchmark record](data/high-depth-read-support-benchmark.json).
These measurements precede the release dependency refresh, which is recorded
separately in [the release catalogue](data/audit-validation-dependencies.json).

## Reproduce

With the `openvax-v1` read bundle already cached and a Topiary checkout containing
its checked-in reference fixture:

```sh
python -m examples.benchmark_high_depth_support /path/to/topiary current.json
PYTHONPATH=/path/to/baseline python /path/to/current/examples/benchmark_high_depth_support.py /path/to/topiary baseline.json
```

Compare the complete `evidence` objects, not only counts. `--variant-id` selects
another checked fixture; `--profile profile.stats` records a separate cProfile
run whose timings include profiling overhead. The benchmark never changes the
read policy or limits discovery to make the comparison faster.
