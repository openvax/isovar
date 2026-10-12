# Index frozen candidate blocks for support matching

Issue [#485](https://github.com/openvax/isovar/issues/485) identified a Cartesian
scan in `_PathSpans`: each read block visited every candidate block on its
contig/strand, even when their genomic intervals were disjoint. Coordinate runs
were already cached; the initial missing-cache diagnosis was incorrect.

The index sorts candidate blocks by genomic start and stores prefix maximum
ends. Two binary searches exclude disjoint blocks before applying the existing
overlap and path-span checks. Repeated, nested and reordered placements remain
eligible. This changes neither sequence matching nor original support identities.
There is no new user parameter.

On a synthetic 4,000-base path split into 400 blocks, 1,000 distinct observations
took 19.96 seconds with the preceding matcher and 0.468 seconds with the index
(42.7× for this component). A seeded 1,000-case comparison also matched the
preceding implementation, including mirrored, repeated, reordered and unplaced
coordinates. Reproduce from this checkout with:

```sh
python -m examples.benchmark_sv_spans --baseline PATH_TO_ISOVAR_1_53_0
```

The complete short-read IGKC benchmark at the unchanged 50,000-record discovery
limit took 77.87 → 71.46 seconds. Complete parsed results, including all array
ordering, were identical. These are single local runs with concurrent validation
work, not a statistically controlled estimate. The 42.7× component improvement
does not imply that speedup for an entire reconstruction. Raw measurements and
engine/input pins are in [the benchmark ledger](span-index-benchmarks.json).

The large ONT recovery uses a separately pinned attempt with unchanged
scientific settings. Its prior four-hour timeout has unknown candidates and
support, so it cannot supply a completed-result runtime or support comparison.
