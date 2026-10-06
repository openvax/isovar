# Complete the original pilot's acquisition-failed SV targets

Follow up [#436](https://github.com/openvax/isovar/issues/436), using the source
date validation from [osteosarc #113](https://github.com/iskandr/osteosarc/pull/113).

1. Derive the exact 22 source/geometry pairs (44 orientations) marked
   unassessable in the frozen 400-view pilot. Retain all original outcomes and
   the historical 13-target recovery as provenance. Reconcile later local
   results and reuse only checksum-verified inputs with the same requested
   regions/source identity; vanished inputs require explicitly fresh acquisition.
2. Use the frozen inventory/transcript-reference snapshots. After upstream
   0.15.5 ships, acquire one source/geometry at a time through its public indexed
   extraction/recovery API. Keep every requested context record, original tag,
   quality and multiplicity. Record seed/partner budgets and acquisition failures;
   timeout/cap/missing source results are unknown, never zero RNA support.
3. Maintain at least 8 GiB free and run each source serially; ONT and 10x can
   proceed independently with separate source metadata. Preserve durable original reads
   and archives; discard only regenerated validation-restore copies after their
   retained archive identities/checksums are verified. Never download whole
   libraries. Retry into separately identified attempts when budgets change.
   Keep staged coverage snapshots under content-derived names so pending and
   completed states both remain inspectable (#456). Preserve the producing
   code/run pins when reusing an earlier verified acquisition.
4. Reconstruct both orientations with the released dense support engine and
   frozen settings. Full-input candidate support and bounded discovery are
   distinct. Preserve ambiguous splice paths, incomplete partner recovery and
   unknown DNA origins. Check exported sequence translation and selected
   full-interval witnesses independently from original BAM records.
5. Publish a pair/view completion ledger, RNA/protein/support rows and readable
   DNA-event/transcript anatomy, including unresolved outcomes. Integrate the
   new selected results with the existing neoORF report and save portable pinned
   bundles. The full-catalogue denominator remains unchanged.
6. Version bump, lint, full tests, green CI, merge and clean-master deployment
   are required. Close #436 only when all its original acquisition-failed pairs
   are accounted for with truthful revised outcomes.
