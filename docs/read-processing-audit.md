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
