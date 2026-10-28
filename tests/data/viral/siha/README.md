# Original SiHa-1 RNA regression slice

These are original public RNA observations from **GSM4785628 / SRR12650977**,
GSE158033 / PRJNA663801, deposited by Zhi-Ming Zheng's NCI laboratory. The
[GEO record](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM4785628)
identifies total-RNA sequencing of the SiHa cervical-cancer cell line.

`original-rna.bam` retains 46 original alignment records from 23 query-name
groups. The earlier experiment used a one-million-pair prefix and competitive
alignment to GDC GRCh38.d1.vd1. Its exact source BAM, full FASTA and replay
hashes are in `manifest.json`; `source-header.sam.gz` preserves the complete
alignment header and minimap2 command/version. No read sequence, quality,
mate flag, tag or selected alignment was changed. Only unused `@SQ` header
entries were removed. Canonical SAM hashes verify every retained record.

Selection uses the saved full-input witnesses: include every query-name group
with at least one segment overlapping the existing 344-nt RNA path and no
recorded placement outside HPV16, then retain **all** original BAM records for
those names. The selection is intentionally conditioned on this one case. It
does not estimate viral abundance or test whole-library specificity/absence.
Three compatible fragments with competing human placements are outside this
slice; their evidence remains in the separate full-input receipt.

The bounded replay yields seven RNA paths from 23 fragments. The selected
344-nt path reproduces the previously reported 98-aa E7 coding prediction and
its L28F difference from NP_041326.1. Its coding bases are tiled by original
RNA; no single segment spans the CDS. RNA strand and transcript completeness
remain unknown. A complete coding prediction is not evidence of protein
expression or an intact full-length viral transcript.

`full-input-replay.json.gz` is the **unchanged historical receipt** for 26 paths
and 71 target fragments, made with ORFdata 0.29.0. It must not be compared as
though these were counts from the 23-fragment slice. The new replay separately
pins the tested ORFdata 0.39.1 catalogue, HPV16 GenBank source and reference
sequence hashes. A changed catalogue requires an explicit fixture review.

The exact HPV16 FASTA already lives at `../hpv16-GDC.fa.gz`. The small replay
uses `custom` reference validation because unused human contigs are absent.
The original BAM header has no M5 checksum, and the replay preserves that
unverified-BAM flag; the manifest separately pins the historical alignment
reference. This is a subset replay, not a new competitive alignment.

From the repository root:

```sh
gzip -cd tests/data/viral/hpv16-GDC.fa.gz > /tmp/siha-reference.fa
isovar viral --bam tests/data/viral/siha/original-rna.bam \
  --genome-fasta /tmp/siha-reference.fa --sample-id SiHa-1 \
  --alignment-reference custom --contig HPV16 \
  --max-discovery-reads 500 --max-paths 32 --output-dir /tmp/siha-replay
```

To reproduce the extraction from the saved original inputs, run
`python examples/retain_siha_fixture.py --bam SOURCE.bam --replay FULL-REPLAY.json
--output-dir NEW-DIRECTORY`. It validates the original BAM checksum, preserves
all selected SAM fields, reruns the existing reconstruction and records hashes.
No download or realignment is performed.

The [NCBI molecular-data policy](https://www.ncbi.nlm.nih.gov/home/about/policies/#data)
places no NCBI restrictions on molecular-data use or redistribution, while
retaining the submitter-rights caveat. This public molecular-data subset is
distributed with source attribution under that policy; it is not relicensed
as project-authored sequence data. The project's code remains Apache-2.0.
