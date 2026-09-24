# Read eligibility across platforms

Isovar decides which alignment records to use with one set of rules, for small
variants (`run_isovar` and the table commands) and for SVs (`isovar sv-rna`), on
reads from any platform. This page lists the rules, explains what each
platform's quality tags do and do not tell you, and records the audit that
unified them.

## Which records are used

| Record | Default | To change |
|---|---|---|
| Unmapped | Skipped | — |
| Failed vendor QC (SAM flag 0x200) | Skipped | — |
| No read name or no sequence | Skipped | — |
| QUAL length differs from the sequence | Skipped | — |
| Marked duplicate (0x400) | Skipped | `--use-duplicate-reads` |
| Secondary alignment (0x100) | Used | `--drop-secondary-alignments` |
| MAPQ below the minimum | Skipped; the minimum is 1 | `--min-mapping-quality N` |
| No base qualities (QUAL `*`) | Used, with base quality unknown | `--require-base-qualities` |
| Rejected by your own predicate | Skipped | `ReadCollector(read_filter=...)` |

The Python equivalents are the `ReadCollector` arguments `use_duplicate_reads`,
`use_secondary_alignments`, `min_mapping_quality` and
`use_reads_without_base_qualities`. `read_filter` receives each original pysam
record and applies in both pipelines. For the tags you can filter on, and an
example, see [read evidence you can filter on](sv-rna.md#read-evidence-you-can-filter-on).

MAPQ is compared as a number, so STAR's 255 for a unique alignment passes the
minimum. Read metadata keep the raw 255, with the normalized MAPQ as unknown
rather than Q255 ([STAR's parameters](https://github.com/alexdobin/STAR/blob/master/source/parametersDefault)).
The stricter splice-inclusion gate of `isovar sv-rna` has its own
[MAPQ 255 option](sv-rna.md#link-an-intronic-start-to-an-observed-transcript-path).

## What quality tags mean

| Tag | Meaning |
|---|---|
| QUAL | Per-base error estimates |
| MAPQ | Confidence in the alignment's placement, not in its bases |
| `rq` (PacBio) | Predicted accuracy of the whole read |
| `qs` (ONT) | Mean Q score of the whole read |
| `np`, `ec` (PacBio) | How much sequencing evidence went into the read |
| `ic`, `is` (Iso-Seq) | Consensus inputs and associations from later processing |
| `OQ` | Original qualities, before recalibration |

Only QUAL says which individual bases are uncertain. Isovar never substitutes
any other tag, or synthetic scores, for missing QUAL. Original qualities in `OQ`
or an earlier file could only be restored after checking sequence, clipping and
orientation, and Isovar does not do this implicitly.
`isovar.read_metadata.record_evidence` returns a record's SAM, ONT and PacBio
tags, whether original qualities are available, and its flags, without combining
them into a weight. SV exports also keep the `RG` and `PG` headers.

ONT split children and duplex consensus reads can share one signal despite
distinct names. SV support reports that lineage separately
([ONT read lineage](ont-read-lineage.md)), and read and fragment counts are not
independent-molecule counts.

Sources: [SAM tags](https://samtools.github.io/hts-specs/SAMtags.pdf),
[SAM flags](https://samtools.github.io/hts-specs/SAMv1.pdf),
[PacBio BAM](https://pacbiofileformats.readthedocs.io/en/13.1/BAM.html),
[CCS filtering](https://ccs.how/faq/reads-bam.html),
[Iso-Seq tags](https://isoseq.how/isoseq-tags.html),
[Dorado tags](https://software-docs.nanoporetech.com/dorado/latest/basecaller/sam_spec/),
[Dorado duplex](https://software-docs.nanoporetech.com/dorado/latest/basecaller/duplex/),
[minimap2](https://github.com/lh3/minimap2/blob/master/minimap2.1).

## Platforms tested

Real Illumina/STAR, ONT/minimap2 and PacBio alignments exercise sequence,
splicing, indels, mates and missing qualities. Synthetic standard-SAM tests
also cover DNBSEQ, IONTORRENT and records that declare no platform. This checks
compatibility with the formats and the policy. It does not show equal
biological sensitivity or calibrated error rates across technologies.

The packaged Sid bundle has 27,077 unique source/SAM-record pairs: 23,109
Illumina, 3,446 ONT and 522 PacBio. All Illumina and ONT records have QUAL, and
515 PacBio records do not. None has `OQ`, none of the PacBio records has
`rq`/`np`/`ec`, and none of the ONT records has `qs`/`dx`/`pi`. These are
selected test records, not estimates of how often tags are missing in a library.

## The 1.21.7 audit

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
