# Original-read SV junction audit

## Question

The broad SV screens count alignments within 500 bp of DNA breakends. Which
candidate links actually resolve an RNA junction, and how many original
fragments support the junction versus a complete observed path?

## Evidence rules

1. Keep source DNA record IDs and screen counts separate. Repeated breakend
   records and two benchmark truth sets may describe the same RNA witness.
2. For an oriented adjacency, run `reconstruct_sv_rna` on the original RNA
   alignments with MAPQ >= 20, actual primary/supplementary records and dense
   support. Report the best event junction's direct fragment count and the
   highest full-path witness count separately. A mate pair is context only.
3. Freeze a small original-read fixture for the three HCC1395 events that
   survive this stricter check. The existing COLO829 selected-site BAM remains
   the paired-read fixture.
4. Keep `references=[]` in this read-geometry audit. It cannot transfer a
   known coding frame or call an SV-derived protein. Ordinary splice
   alternatives and DNA-to-RNA causation remain unresolved.

## Deliverables

- A reproducible script and TSV of audited candidate events, with explicit
  BAM identity and limits.
- A regression test against source-derived original reads, including a
  low-MAPQ COLO829 near-join that must not be called at MAPQ >= 20.
- Report text and PDFs showing direct-junction and whole-path counts with
  their distinct meanings.
- [Isovar #553](https://github.com/openvax/isovar/issues/553) for an
  intermediate, path-specific local-context metric;
  entire-path exact support is too strict for noisy long reads, while shared
  junction support is too broad to assign a transcript isoform.
