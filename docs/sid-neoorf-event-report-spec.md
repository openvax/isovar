# DNA → RNA → ORF report

Extend PR #447's focused candidate screen with an event-by-event report covering
all nine selected SV geometries and all 32 indel/splice nominations, including
the unresolved MUC3A allele and zero-candidate outcomes. Preserve distinct
geometries, transcript interpretations, RNA products and timepoints.

For each event, retain original DNA call/allele identifiers, caller/sample,
coordinates, REF/ALT and available DNA support. RNA-only nominations must say
that the DNA event is unknown. Describe retained transcript regions, exon/
intron/UTR/CDS placement, RNA joins and coding-frame consequences using the
pinned Ensembl 115 models. Unresolved events still get transcript anatomy and
event locations, without inventing a mature transcript or translated protein.

Provide every screened SV ORF and frameshift coding-window sequence, with stop
status, transcript IDs, exact sequence/window support by product, quality and
barcode limitations, and source/result pins. Keep in-frame coding controls in
the complete appendix, explicitly labeled. Distinguish DNA caller support,
RNA allele/junction support, complete-window witnesses and gene expression.
Do not sum overlapping hypotheses, geometry aliases or processing products.

Publish a readable event report plus complete sequence/support appendices and
an exportable ledger. Verify counts and references against the frozen screen;
no new whole-library acquisition. Re-run required local gates and final-head
CI, then merge/deploy once CI is green.

Validation adjustment (#450): the full run exposed two mocked-acquisition tests
depending on actual host free disk space. Control that value in the offline
fixture and test low-space rejection/the exact 8 GiB boundary explicitly. Keep
the production acquisition guard unchanged; re-run the final gates.

CI adjustment (#451): all three Python jobs exposed the original-BAM checker's
undeclared SAMtools CLI prerequisite. Install SAMtools in CI and document the
requirement; keep the native verification/test enabled and repeat final-head CI.
