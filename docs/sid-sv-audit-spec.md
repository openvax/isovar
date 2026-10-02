# Sid catalogue-wide SV RNA reconstruction

## Objective

Run Isovar's existing SV RNA reconstruction across Sid's nominated structural
events and all relevant RNA alignment products. Produce a reproducible ledger
in which every intended event and source has an outcome, including unavailable
data, unresolved geometry, no observed junction, ambiguous paths and bounded
searches. Rank illustrative candidates by actual ORF support, separately from
junction support. This does not establish translation or clinical relevance.

## Inputs and identity

- Freeze the shared osteosarc SV candidate catalogue, including single breakends
  and aliases, and the broader fusion survey. Preserve original calls and source
  checksums; use an explicit coordinate conversion for each source format.
- Freeze the RNA source inventory, including short-read, ONT, PacBio and 10x
  products across available timepoints. Record reference assembly, assay,
  library/specimen assertions and processing provenance. Unknown identities stay
  unknown. Reprocessed BAMs are separate observations of a library, never added
  as independent support.
- Preserve existing regression fixtures. Acquire indexed regions and linked
  original records through osteosarc, with receipts, hashes and explicit limits.
  Never scan or download a whole remote BAM as a fallback.

## Reconstruction and output

- Reuse `reconstruct_sv_rna` with both justified orientations, pinned transcript
  references and recorded parameters. Keep alternative paths and proteins.
- Separate local junction evidence, assembled sequence support, and direct
  full-ORF witnesses. Retain native-query orientation, missing quality, read
  groups, cell/UMI labels, complete/partial ORFs and initiation uncertainty.
- Do not silently skip unresolved geometry, reference mismatches, missing
  indexes, empty inputs, failed acquisitions or exhausted reconstruction caps.
- Cache immutable inputs and resumable per-event/source results. Verify resumed
  inputs and parameters. Provide an offline rerun and a coverage ledger that
  fails if any intended combination is missing or duplicated.

## Validation and delivery

Independently check translation and original-read support in representative
positives and controls. Cover catalogue normalization, complete accounting,
aliases/reprocessed products, missing qualities, unsupported events and failed
or truncated work. Run the actual catalogue audit and publish its precise
denominators, results and limitations. Lint, full tests and GitHub CI must pass;
land via a versioned PR and deploy using the repository release workflow.

Primary input documentation: [shared SV catalogue](https://iskandr.github.io/osteosarc/sv-candidates/),
[fusion survey](https://osteosarc.com/fusions/), and
[SAM specification](https://samtools.github.io/hts-specs/SAMv1.pdf).

## Scale check and revised execution scope

The first 50-geometry PacBio batch fetched 404,697 original context records.
Partner recovery reached its wall-time bound. Reconstructing unrelated regional
junctions also consumed substantial time without adding event-crossing ORFs.
Before expanding to all products, benchmark an explicit event-focused seeding
option. It must keep the same collected context, overlap extension and ambiguous
event-compatible splices, and report excluded unrelated seeds. The general API
default remains the broader search. Acquisition failures stay visible; do not
turn timeout-truncated inputs into negative evidence.

The same 50-geometry batch exceeded the 500,000-record acquisition cap in ONT
T2, bulk RNA T2 and 10x T2. Partition oversized batches recursively, preserving
the original windows, transcript context and intended geometry matrix. Store
the partition tree and bind each result to its leaf receipt. Only a single
geometry still exceeding the cap receives a terminal record-limit failure.
An interrupted partition stays pending, and missing children cannot pass the
coverage check. This changes request scheduling, not which reads qualify.

Publish each immutable checkpoint atomically from a temporary file in the same
directory, without replacing different existing content. Apply the per-view
deadline only to reconstruction and summarization; persist after cancelling the
alarm. Interrupted writes must leave no final checkpoint, and timeouts must not
retain a partially computed candidate summary (review issue #435).

The four-platform pilot also exposed duplicated upstream query-name lists in
acquisition receipts and command input files (osteosarc issue #107). Measure and
deduplicate this metadata before allocating storage for the full catalogue;
preserve original read hashes, semantic request identities and all provenance.
