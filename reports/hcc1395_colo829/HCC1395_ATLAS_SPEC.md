# HCC1395 whole-run RNA atlas: implementation plan

Issue: [Isovar #533](https://github.com/openvax/isovar/issues/533).

1. Pin the SEQC2 v1.2.1 matched HCC1395/HCC1395BL call set, the
   GRCh38 no-alt reference, and the complete original ONT direct RNA run
   SRR23456709. Validate the archive and retain indexed original-read
   alignments at every SEQC2 SNV/indel and both ends of each Talsania SV.
2. Join every source small-variant record to reference, alternate, and other
   RNA read counts. Report no coverage separately from observed reference.
   Evaluate every alternate-supported call for a local Isovar protein context,
   with explicit no-translation outcomes. Keep matched-normal DNA evidence
   distinct from RNA in a separate culture.
3. Screen the Talsania SV source records for breakend-window RNA coverage and
   same-read split-alignment candidates. Preserve record and fusion-source
   identity; do not infer a junction or protein from proximity alone.
4. Commit bounded original-read fixtures for offline positive, reference-only,
   no-coverage, and SV screen tests. Publish machine-readable tables and a
   versioned narrative with exact tool versions, checksums, and limitations.

The direct RNA assay has higher error rates than the COLO829 Illumina assay.
An ALT observation is an allele witness, not a validated somatic transcript.
Local protein contexts require separate read and reference-transcript
compatibility; none asserts a full-length isoform.
