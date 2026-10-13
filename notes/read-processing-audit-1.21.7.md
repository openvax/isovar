# Read-processing audit for 1.21.7

Maintainer working note, split from the published guide [Read eligibility across platforms](../docs/read-processing-audit.md).

Before 1.21.7, small-variant collection accepted records failing vendor QC while
SV collection rejected them. The audit set out to:

- apply one eligibility policy to both pipelines, covering vendor QC, missing
  QUAL, duplicate and secondary flags, and MAPQ;
- keep unknown qualities unknown, and keep native evidence tags separate;
- extract ONT, PacBio and general SAM evidence through one shared helper, without
  parsing optional tags on every read during ordinary collection;
- offer a caller-supplied record filter to both pipelines;
- stop allele-only logging from reading gene annotation
  ([#295](https://github.com/openvax/isovar/issues/295));
- measure collection on the existing osteosarc-derived fixtures, adding no read
  data.

All of these shipped in 1.21.7. The mate merger also stopped sorting singleton
groups and scanning all their bases for a sort key. Tag extraction stays outside
ordinary collection and does not materialize large auxiliary arrays.

### Measurements

The [benchmark script](../examples/benchmark_read_processing.py) decodes only the
packaged osteosarc-derived records. It runs 51 locus queries: 32 Illumina, 16 ONT
and 3 PacBio. Records can recur across queries, so this is not a depth or
support estimate.

Median CPU time over 21 alternating runs per case and mode, 1.21.6 (commit
`13cfbd778a5b14a2f63846d2c29a103d8cc918a3`) → 1.21.7:

| Platform | Queries | Public collection | Compact collection |
|---|---:|---:|---:|
| Illumina | 32 | 102.78 → 101.96 ms | 99.46 → 99.87 ms |
| ONT | 16 | 112.04 → 101.15 ms | 97.88 → 97.63 ms |
| PacBio | 3 | 8.58 → 7.56 ms | 7.39 → 7.39 ms |

- All six observation fingerprints match, retaining 1,140, 784 and 70
  observations respectively. Peak allocation during collection is essentially
  unchanged.
- Ordinary collection stays level. The public long-read path gains from no
  longer scanning singleton sort keys.
- The optional, expanded metadata extraction takes 11.40 ms for 2,061 Illumina
  records, 11.49 ms for 979 ONT records and 7.70 ms for 781 PacBio records. It
  exports more fields than the earlier, shorter list and costs more.

These are local microbenchmarks. They exclude acquisition, SAM decoding,
annotation and protein reconstruction. Separate-process wall times fluctuated
too much to compare, so runs alternate between implementations and record both
wall and CPU samples. CI makes no timing assertions. Full samples, allocation
peaks, input and script fingerprints, and Python and platform details are in
[read-processing-benchmark.json](read-processing-benchmark.json). To rerun:

```sh
python -m examples.benchmark_read_processing output.json \
  --baseline /path/to/clean/isovar-1.21.6-checkout --repeats 21
```
