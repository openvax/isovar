# Original Sid long-read records at nominated fusion breakpoints

Unchanged SAM records from the public Sid Sijbrandij osteosarcoma bucket,
retrieved 2026-09-18 by `build_long_read.py`. The script queries each
event's breakpoint windows (±1000 bp) and every supplied reference exon with
indexed `samtools view -M`. It keeps segments (read group, QNAME, mate) that
have records overlapping **both** breakpoint windows: split and spliced-through
junction reads. Reads with only one partner aligned are not included. Breakpoints and
Ensembl 87 models come from the adjacent `corpus`/`coding-corpus` inputs. The
manifest pins source URLs, record counts and file hashes.

| File | Source | Exercises |
|---|---|---|
| `TPST1--CRCP.PacBio-T1` | T1 Kinnex/Iso-Seq, `pbmm2 --preset ISOSEQ` | pbmm2 places the junction's 8-nt homology on CRCP (assignment [3, −3]): the same adjacency as ONT's |
| `FOXO3--STRADA-CCDC47.ONT-T1-tagged` | T1 ONT single-cell, minimap2 | noisy reads support their own junction without matching the consensus end to end |
| `ATP5MG--KMT2A.PacBio-T1` | T1 Kinnex/Iso-Seq | exact join between an annotated ATP5MG exon end and KMT2A exon start: read-through-ambiguous |

Dataset facts relevant to reuse:

- PacBio T1 and ONT T1 share ~75–79% of cell barcodes (and none with ONT T2):
  one single-cell cDNA library sequenced on two platforms, not independent samples.
- The ONT `*_dedup.bam` products keep one alignment record per molecule for
  TPST1–CRCP split reads, removing the junction: 121 T1 molecules in `.tagged.bam`, 0 in
  `_dedup.bam`. Use `.tagged.bam` for split-read fusions; count molecules by `CB`/`UB`.
- The ATP5MG–KMT2A join is also what ordinary read-through splicing of adjacent,
  same-strand genes produces. RNA alone does not establish a rearrangement.

Sources: [SAM](https://samtools.github.io/hts-specs/SAMv1.pdf),
[SA tag](https://samtools.github.io/hts-specs/SAMtags.pdf) and the
[osteosarc fusion catalogue](https://osteosarc.com/fusions/).
