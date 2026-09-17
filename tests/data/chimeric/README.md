# Original osteosarc split-read records

Two observed paths from the public Sid Sijbrandij RNA data, retrieved 2026-09-16.
These exercise phasing provenance, not coding-fusion reconstruction. The tests
probe synthetic alleles at observed bases; they do not assert called mutations.
SAM records (including CIGAR, flags, sequence, qualities and tags) are unchanged.
Headers retain only the referenced chromosomes and read groups; sort order is
unspecified. No supplementary record is synthesized from an SA tag.

## ONT RNA

- Fixture: `osteosarc-ont.sam`
- Source: https://sid-sijbrandij-osteosarc-dataset.s3.us-west-2.amazonaws.com/ONT/IPISRC044_ONT_upload/IPISRC044_ONT/processed/IPISRC044_T1_sclrs_ONT/IPISRC044_T1_sclrs_ONT/IPISRC044_T1_sclrs_ONT_dedup/IPISRC044_T1_sclrs_ONT_dedup.bam
- Acquired regional BAM SHA-256: `71be5e3681c02e86155df8ae25cb0f91814e3d94a741a50476b3d3e6c19bb923`
- QNAME: `46de198c-6f57-490f-9f0e-1679fb487660_0`

## Short-read RNA

- Fixture: `osteosarc-short.sam`
- Source: https://sid-sijbrandij-osteosarc-dataset.s3.us-west-2.amazonaws.com/kamil/oncoanalyser/IPISRC044_T1_ucla/alignments/rna/IPISRC044_tumor_T1_ucla_rna.md.bam
- Acquired regional BAM SHA-256: `dc5e73d0004f7402b76602f6d1135876fc772b8c36675af68410d7baa0910272`
- QNAME: `A01231:888:HGH32DSXC:1:1535:20898:9909`

The ONT path crosses chr6/chr17 and includes a reverse-strand hard-clipped
supplementary record with `356H5M1I128M1H`. Its partner's SA tag reports the
approximation `356S133M1I1S`, not an exact base mapping. The short-read fixture
also retains the distinct paired mate (without SA); the split-path regression
selects only the two linked records of the other segment.

Path semantics follow [SAM](https://samtools.github.io/hts-specs/SAMv1.pdf)
and [SA tags](https://samtools.github.io/hts-specs/SAMtags.pdf).
Minimap2's [SA writer](https://github.com/lh3/minimap2/blob/master/format.c)
uses query/reference spans to summarize CIGARs. Never use that approximation
to extract junction bases or construct a missing alignment.
