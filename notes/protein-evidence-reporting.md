# Protein evidence semantics (#590, #551)

Preserve contributor counts and marked mutation intervals as existing API
metrics. They do not establish exact coding coverage or changed amino acids.

1. Extend the shared coding-interval assessment with an explicit merged-mate
   span subset and its witness IDs. Keep exact coding interval, per-base depth,
   single-observation spans and observed-stop coverage together. Do not pool
   different nucleotide paths or claim observed translation.
2. Centralize sequence comparisons in the library. Keep annotation-reference
   comparison separate from the same-RNA focal-allele counterfactual. Export
   identical/divergent/unresolved states, verified differences and marked span
   independently. Preserve alternate models and provenance. A synonymous focal
   allele can coexist with other unexplained nonreference residues.
3. Make the generic report consume protein-hypothesis exports, and use the same
   classification/selection helpers in cohort reports. Separate context counts
   from changed-protein counts and show uncertainty rather than converting old
   tables without comparison provenance to verified changes.
4. Repair the historical focused report's contributor-as-complete-window prose
   and FASTA headers; use shared assessments for a pinned original-RNA recount.
   Cover CDKN2A, FHL2, GOLGA5 and the C8orf76 control. Re-render and visually
   inspect figures/PDFs; retain original source and interval identities.

Validation: synonymous ALT/reference-identical sequence, an independently
changed RNA background, a focal deletion gap, frameshift and stop boundaries,
partial contributors, merged mates, source identity deduplication, complete
versus local coding candidates, and ambiguous or missing comparators. Run the
complete suite and CI before merge and clean-master deployment.

Reference semantics: SAM v1.6 distinguishes templates and sequenced segments
(https://samtools.github.io/hts-specs/SAMv1.pdf); translation uses the recorded
NCBI genetic code (https://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi).
Neither exact RNA coverage nor a reference frame establishes translation.
