# Index frozen candidate blocks for support matching

## Question, result, and code decision

**Question:** Did repeated scans of disjoint candidate blocks dominate exact
SV support matching, and could those scans be avoided without changing
original-read witnesses?

**Result:** A synthetic 4,000-base, 400-block component test improved from
19.96 to 0.468 seconds (42.7×). The complete short-read IGKC replay changed
from 77.87 to 71.46 seconds with identical parsed results. Subsequent witness
and shared-query refinements preserved those results but did not establish a
further complete-run speedup. These are historical single-machine measurements,
not a current performance guarantee.

**Code influenced:** [PR #484](https://github.com/openvax/isovar/pull/484)
added the `_PathSpans` interval index to exclude disjoint blocks;
later changes reject incompatible witness offsets early and share exact
interval searches and equal witness rows within bounded batches. The code
retains every original read identity, placement and repeated occurrence.
There is no new user parameter. See [issue #485](https://github.com/openvax/isovar/issues/485).

## Methods and successive refinements

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
engine/input pins are in [the benchmark ledger](sv-span-index-485.json).

The large ONT recovery uses a separately pinned attempt with unchanged
scientific settings. Its prior four-hour timeout has unknown candidates and
support, so it cannot supply a completed-result runtime or support comparison.

## Exact witness-query refinement

The subsequent refinement rejects incompatible read/path query offsets early,
indexes the disjoint union of candidate genomic blocks for the initial CIGAR
overlap screen, and compares overlapping linear placement runs when checking
witnesses. Sequence searches and repeated interval checks use bounded internal
caches. Cache exhaustion retains every occurrence through the original exact
calculation. Both strands, unplaced bases, repeated placements and original
witness identities keep their existing semantics.

The complete short-read benchmark remained identical after this refinement:
71.46 → 72.59 seconds and approximately 1.48 GiB peak RSS in both runs. These
single shared-machine measurements do not establish a short-read speedup.
The [witness-query ledger](sv-witness-query-485.json) records full inputs,
engine pins and result equality. Validation includes 5,919 full tests and 206
differential tests against the preceding witness and complete-support code.

## Share interval searches and witness rows

Within each bounded support batch, identical candidate DNA intervals now share
exact searches. Search only direct-junction observations; retain every original
identity and every repeated occurrence. Reject incompatible placements within
the requested interval before scanning for flank anchors. Share complete equal
witness rows in memory, including observation identity and all junction links.
Internal cache limits fall back to exact computation without limiting evidence.

In one 300-template ONT component batch, complete retained state matched the
preceding implementation. Three alternating comparisons gave a 2.74× median
speedup; retained Python objects occupied 4.98 → 1.80 MB. All 2,326 witness rows
remained present, backed by 460 distinct equal row objects. These are component
measurements, not a complete-ONT speedup or process-RSS claim. The complete
short-read result remained identical (77.70 seconds in this shared-machine run;
no short-read speedup claim). See the [comparison ledger](sv-shared-query-485.json).
