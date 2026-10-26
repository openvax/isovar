# COLO829 SV3: regional RNA evidence for a reading-frame anchor

The result is **no event-linked frame anchor in the scanned inputs**.
The original-read search found exon-linked sequences in the same
downstream reference block as the SV3 RNA join, but none of those
spots contains the SV3-specific 43-base junction. A shared reference
segment does not phase two upstream RNA branches into one transcript.
This advances [#563](https://github.com/openvax/isovar/issues/563).
The exploratory ORFs remain separate from assigned proteins.

The [machine-readable audit](sv3_frame_anchor_audit.json) scanned
**52,212 alignments** from
**26,976 SRA spots** in the [pinned CNIH3
candidate BAM](../../tests/data/colo829_whole_run/SRR8615617.cnih3_source_candidates.grch38.bam),
without a MAPQ threshold. Exact matches require at
least 20 path bases; exon alignments require 15 aligned bases. The
ten raw ends of the five known SV3 junction spots were also checked.
Across both sources, 1,242 spots
have at least 20 exact bases of the selected joined path; most matches
are to ordinary reference sequence away from the SV3 junction.
The BAM contains 6,233
spots with a qualifying Ensembl 87 exon alignment, including
673 with a coding-exon
alignment. Coding-region expression alone does not link an RNA
molecule to the SV3 join.

## Exon-linked candidates

**9 spots** have both a versioned
CNIH3 exon alignment and at least 20 exact bases of the selected SV3
path; **6** reach 30 bases.
Every matched segment is in the ordinary downstream reference block,
starting 280 bases after the 19 unplaced SV3 join bases.
These are
candidate alternative exon-to-region RNA links, not reads spanning
the SV3-specific join.

| SRA spot | Best exact path match, GRCh38 1-based | Exon end(s), 1-based | Link | Clip quality | SV3 43-mer |
|---|---:|---|---|---|---|
| `SRR8615617.103128066` | 32 nt at chr1:224598612-224598643 (path `1387:1419`) | 224583267 | same-end clipped | right-clip min Q 29 | no |
| `SRR8615617.104275726` | 30 nt at chr1:224598612-224598641 (path `1387:1417`) | 224583267 | same-end clipped | right-clip min Q 8 | no |
| `SRR8615617.110017311` | 30 nt at chr1:224598612-224598641 (path `1387:1417`) | 224583267 | same-end clipped | right-clip min Q 30 | no |
| `SRR8615617.7298052` | 20 nt at chr1:224598612-224598631 (path `1387:1407`) | 224566264 | same-end clipped | right-clip min Q 30 | no |
| `SRR8615617.7300160` | 22 nt at chr1:224598612-224598633 (path `1387:1409`) | 224566264 | same-end clipped | right-clip min Q 10 | no |
| `SRR8615617.7300195` | 32 nt at chr1:224598612-224598643 (path `1387:1419`) | 224583267 | same-end clipped | right-clip min Q 31 | no |
| `SRR8615617.7303961` | 20 nt at chr1:224598612-224598631 (path `1387:1407`) | 224583267 | same-end clipped | right-clip min Q 31 | no |
| `SRR8615617.7303964` | 31 nt at chr1:224598612-224598642 (path `1387:1418`) | 224583267 | same-end clipped | right-clip min Q 31 | no |
| `SRR8615617.7306533` | 31 nt at chr1:224598612-224598642 (path `1387:1418`) | 224587420 | same-end clipped | right-clip min Q 31 | no |

The exon coordinates above are **exon ends**, not the start of
the downstream clipped segment. The JSON retains each original
alignment, mate, CIGAR, MAPQ, exon model, clip quality, and exact
path interval. Multiple alignments and overlapping transcript exons
do not add independent spot support.
**0** exon-linked spots align to
a coding block in Ensembl 87. A soft-clipped sequence match
does not by itself establish a mature splice or reading frame.

## Event-specific and frame evidence

The exact SV3 43-mer has **5**
distinct SRA spots; **3** have
placed split-junction support. **0**
have any qualifying exon alignment in the BAM, and none is among
the exon-linked path matches.
1 exact 43-mer
spot requires the pinned raw FASTQ end because its sequence is not
an exact 43-mer in the candidate BAM record.
The retained joined assembly has **0**
exact versioned exon-anchor bases and
**0** event-linked
frame translations. Ensembl 87 coding CNIH3 start: `ENST00000272133` at chr1:224,617,175.

| Exploratory ORF | Original spots matching its ATG | Single ends covering ATG to join core / complete 43-mer | Full ORF overlap chain | Annotated or observed initiation |
|---|---:|---:|---:|---|
| 37 aa | 47 | 1 / 0 | 2 spots | no |
| 162 aa | 51 | 2 / 0 | 7 spots | no |

An ATG sequence match is not evidence of translation initiation.
An overlap chain combines different spots and does not establish
one full-length molecule or a selected reading frame. The exact
start-to-join spot IDs are in the [path-fusion audit](sv3_path_fusion.json).
The [polished-isoform cross-check](colo829_sv3_polished_check.md)
found no exact SV3 junction among
68,553 representatives; its
reverse-motif representatives are
`PB.1153.13`, `PB.1153.23`.
counts are a separate assay, not original-read support for this join.

## Scope and next evidence

This BAM retains candidates from an alignment-based screen of the
full SRR8615617 source. Completely unaligned connectors can be missed.
The Isovar path search reached its 80-path limit.
Matched-normal RNA
and original COLO829 PacBio reads were not assessed. The DNA call and
RNA source also came from different COLO829 cultures or passages.
A supported path must connect the exact 19-base SV3 join to a
versioned exon or other established ORF anchor, with orientation and
frame evidence. Overlapping reads can establish that path; one
full-length read is not required. The alternative exon-to-downstream
links here are specific candidates for comparison, not a substituted
SV3 proteoform.

Regenerate offline with
`python -m reports.hcc1395_colo829.build_sv3_frame_anchor_audit`.
