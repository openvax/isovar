# Cross-platform read-processing audit

## Scope and acceptance criteria

- Use one record eligibility policy for small variants and SVs, including
  vendor-QC failure, missing QUAL, duplicate/secondary flags and MAPQ.
- Preserve unknown qualities. No alignment, pass-count, barcode or consensus
  tag can reconstruct per-base QUAL. Keep native evidence separately.
- Extract useful ONT, PacBio and general SAM evidence through a shared helper;
  make a caller-supplied record filter available to both collection paths.
  Default collection must not parse optional metadata on every read.
- Remove annotation access from allele-only logging (#295). Measure collection
  and metadata extraction on the existing osteosarc-derived minimal fixtures;
  add no new read data. Verify actual Illumina, ONT and PacBio reads, plus
  standard SAM records from unknown/other platforms.
- Run regression tests, lint, the full suite and GitHub CI before release.

## Scientific definitions

QUAL contains per-base error estimates; MAPQ describes alignment placement.
`rq` is predicted read accuracy, while ONT `qs` is a mean read Q score.
`np`/`ec` count sequencing evidence; Iso-Seq `ic`/`is` describe later consensus
inputs/associations. None specifies which individual bases are uncertain.
Original qualities may exist in `OQ` or an earlier sequencing file, but recovery
requires checking sequence, clipping and orientation. Isovar does not replace
QUAL with OQ or synthetic scores implicitly.

Sources: [SAM tags](https://samtools.github.io/hts-specs/SAMtags.pdf),
[SAM flags](https://samtools.github.io/hts-specs/SAMv1.pdf),
[PacBio BAM](https://pacbiofileformats.readthedocs.io/en/13.1/BAM.html),
[CCS filtering](https://ccs.how/faq/reads-bam.html),
[Iso-Seq tags](https://isoseq.how/isoseq-tags.html),
[Dorado tags](https://software-docs.nanoporetech.com/dorado/latest/basecaller/sam_spec/),
[Dorado duplex](https://software-docs.nanoporetech.com/dorado/latest/basecaller/duplex/),
[minimap2](https://github.com/lh3/minimap2/blob/master/minimap2.1).

## Findings and changes

Isovar 1.21.6 accepted SAM flag 0x200 in small-variant collection while SV
collection rejected it. The shared eligibility method now rejects vendor-QC
failures in both pipelines. Both accept missing QUAL by default and keep the
existing duplicate, secondary and numeric MAPQ options. STAR's default MAPQ
255 remains accepted; exported metadata reports its raw value and an unknown
normalized MAPQ. This preserves unique STAR alignments without calling them
Q255. See [STAR's parameters](https://github.com/alexdobin/STAR/blob/master/source/parametersDefault).

`ReadCollector(read_filter=...)` adds an optional original-record predicate to
both pipelines. `isovar.read_metadata.record_evidence` exposes native SAM,
ONT and PacBio tags, original-quality availability and flags without assigning
a combined weight. SV exports also retain RG/PG headers. See
[the field definitions and filtering example](sv-rna.md#read-evidence-you-can-filter-on).

Allele collection no longer accesses gene annotation solely for logging (#295).
The mate merger avoids sorting singleton groups and scanning all their bases
for a sort key. Selected tag extraction stays outside ordinary collection and
does not materialize large auxiliary arrays.

The existing packaged Sid bundle contains 27,077 unique source/SAM-record pairs:
23,109 Illumina, 3,446 ONT and 522 PacBio. All Illumina and ONT records have QUAL;
515 PacBio records do not. None has OQ; none of the PacBio records has rq/np/ec,
and none of these ONT records has qs/dx/pi. These are selected test records,
not library-wide missingness estimates. No read fixtures were added or enlarged.

Real Illumina/STAR, ONT/minimap2 and PacBio alignments exercise sequence,
splice, indel, mate and missing-quality behavior. Synthetic standard-SAM tests
also cover DNBSEQ, IONTORRENT and an absent platform declaration. This verifies
format/policy compatibility, not equal biological sensitivity or calibrated
error rates across technologies.

ONT split children and a duplex consensus can share original signal lineage
despite distinct read names. SV support reports that lineage separately
([ONT read lineage](ont-read-lineage.md)); segment/fragment counts are not
independent-molecule counts.

## Measurements

The [reproduction script](../examples/benchmark_read_processing.py) decodes only
the packaged osteosarc-derived records. It runs 51 locus queries: 32 Illumina,
16 ONT and 3 PacBio. Input records may recur across queries;
this is not a new sequencing-depth or support-count estimate.

Compared with 1.21.6 (commit `13cfbd778a5b14a2f63846d2c29a103d8cc918a3`),
21 alternating measurements per case/mode produced these median CPU times:

| Platform | Queries | Public collection, before → after | Compact collection, before → after |
|---|---:|---:|---:|
| Illumina | 32 | 102.78 → 101.96 ms | 99.46 → 99.87 ms |
| ONT | 16 | 112.04 → 101.15 ms | 97.88 → 97.63 ms |
| PacBio | 3 | 8.58 → 7.56 ms | 7.39 → 7.39 ms |

All six complete observation fingerprints match; retained counts are 1,140,
784 and 70, respectively. Collection allocation peaks are effectively unchanged.
The optional, expanded metadata extraction costs 11.40 ms for 2,061 Illumina
records, 11.49 ms for 979 ONT records and 7.70 ms for 781 PacBio records. It grows
the exported metadata and costs more than the former, smaller field list.

These are local microbenchmarks, excluding acquisition, SAM decoding, annotation
and protein reconstruction. Ordinary collection stays essentially level; the
public long-read path benefits from avoiding singleton sort-key scans. Earlier
separate-process wall timings fluctuated substantially, so the reported
comparison alternates implementations and records both wall and CPU samples.
There are no machine-dependent timing assertions in CI. Full samples, allocation
peaks, input/script fingerprints and Python/platform details are in
[read-processing-benchmark.json](read-processing-benchmark.json).

```sh
python -m examples.benchmark_read_processing output.json \
  --baseline /path/to/clean/isovar-1.21.6-checkout --repeats 21
```
