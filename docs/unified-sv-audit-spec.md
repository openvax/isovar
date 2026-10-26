# Unified nominated-SV RNA audit (#305, #541, #563)

## Goal

Use the existing `reconstruct_sv_rna` engine for a reproducible, source-separated
audit from nominated DNA events through alternative RNA paths to justified
protein candidates. SV3 is an unresolved-frame regression case; HCC1395
EIF3K–CYP39A1 is the positive coding-fusion control. A completed computation
must not imply a biologically resolved event.

## Implementation

1. Assess every retained path in a library API: original junction witnesses,
   local branch support, exon/CDS links, competing placements, exact sequence
   coverage, frame evidence, and partial/complete proteins remain separate.
   Move reusable SV3 report calculations into this API; report scripts render
   its results and retain pinned source-specific provenance.
2. Reconcile nucleotide paths with source-identified Varcode predictions,
   ordinary transcripts and supplied alternative placements. Preserve support,
   informative contradiction, indistinguishability, insufficient coverage and
   unassessed outcomes. Sequence overlap and mate locations do not establish
   phase or a reading frame. Reuse existing frame and supplied-fusion APIs.
3. A versioned manifest pins event identity, oriented inputs, reference models,
   source/library/culture, file hashes and discovery budgets. The same runner
   audits the separate COLO829 SV inventories, original tumor short/long reads
   and available normal RNA. Missing data and limited searches are explicit;
   no-candidate rows remain in the cohort report. Keep small-variant atlas
   results linked rather than recomputing them through the SV engine.

## Verification and delivery

Test all-path assessment, evidence discrimination, source identity, repeated
alignments, reverse orientation, partial/full coding paths, unresolved SV3 and
the original-read positive fusion. Reproduce committed cohort results from
bounded fixtures with pinned hashes. Run lint, the full suite and GitHub CI;
merge and deploy from a clean master. Update #541 to reflect its completed
atlas, consolidate overlapping requirements under #305, and close issues only
when their stated acceptance criteria are satisfied.

## Scientific basis

SAM v1.6 distinguishes templates, segments and alternative/supplementary
alignments: https://samtools.github.io/hts-specs/SAMv1.pdf . Competitive fusion
assessment should retain ordinary-transcript explanations and compare tumor
and normal evidence: https://pmc.ncbi.nlm.nih.gov/articles/PMC10261907/ .
