# Sid read-end regression records

`osteosarc.sam` contains four unchanged primary SAM records selected from the
cached T1 BAM subsets below. The header is reduced to their reference sequences
and original read groups. Sequence, flags, CIGAR, qualities (including missing
QUAL), tags and read names are unchanged. These are not simulated adapter reads.

- ONT: `a81a7468-ebc9-440f-a8de-e5f6c6591852_0` and
  `231d999d-c208-47d7-9133-57b07ae7f5cf_0`. Their terminal clips are `CT` + 24 A's
  and 26 A's. The latter has five additional A bases inside the alignment: those
  must not be trimmed. Chemistry and biological full-length status are unknown.
  [Source BAM](https://sid-sijbrandij-osteosarc-dataset.s3.us-west-2.amazonaws.com/ONT/IPISRC044_ONT_upload/IPISRC044_ONT/processed/IPISRC044_T1_sclrs_ONT/IPISRC044_T1_sclrs_ONT/IPISRC044_T1_sclrs_ONT.tagged.bam)
- Illumina read 2: `A01231:888:HGH32DSXC:3:1373:1163:22889`, `5S105M41S`.
  The right clip contains the reported Illumina adapter motif; matching the
  TruSeq R2 trimming sequence requires one edit. This is a sequence-consistent
  adapter assignment, not proof of the preparation's kit identity.
  [Source BAM](https://sid-sijbrandij-osteosarc-dataset.s3.us-west-2.amazonaws.com/kamil/oncoanalyser/IPISRC044_T1_ucla/alignments/rna/IPISRC044_tumor_T1_ucla_rna.md.bam)
- PacBio `molecule/10052202`, RG `e4927d21`: a processed group-deduplicated
  transcript with absent QUAL. Retained as a negative/preservation test, not a
  validated adapter-positive or full-length transcript.
  [Source BAM](https://sid-sijbrandij-osteosarc-dataset.s3.us-west-2.amazonaws.com/pacbio/IPISRC044_T1_sclrs_live_pbmm2_mapped.bam)

Acquisition receipts, unchanged subset BAMs and SHA-256 manifests are retained in
the local `figures/osteosarc/2026-09-18_13-26-15-937302Z/sv-realignment/` evidence
bundle. See [#302](https://github.com/openvax/isovar/issues/302) for the investigation
and the broader platform-aware design. Regression tests require no network.
