# Shared cohort RNA analysis and reports

One `isovar cohort --manifest ... --output-dir ...` command runs source-separated
small-variant/ncORF reconstruction and the existing SV, supplied-fusion and viral
APIs, then renders per-sample reports from their results. The manifest names all
inputs, references, library identities, unavailable sources and analysis budgets.
An inventory is not a claim that every input has been assessed successfully.

Small variants retain every call, original allele support, all retained protein
and nucleotide alternatives, exact coding-interval coverage, validated ncORF
competitors and explicit errors. The supplied-ORF path is shared with the library.
Contributors, merged-mate full spans, single-original-read spans and per-base
tiling remain separate (#590). Reference sequence comparisons do not fill RNA.

Resume is input-, parameter- and implementation-pinned, with per-call checkpoints;
failures retain their stage and message and cannot become zero-support rows.
Libraries, timepoints and cultures remain separate. Reprocessed products must be
declared and are never pooled automatically. All original data provenance is
retained, including exclusions and unknown sample assignments.

Reports show source coverage and missing analyses first, then variant landscapes,
support funnels and top examples by category. Sequence panels compare observed
protein with reference, use distinct focal/germline/co-somatic/unclassified
colors, mark deleted residues and show coding coverage. Figures and full data
exports accompany the PDF. PDF rendering can be rerun without reconstruction.

Validation includes original COLO829 merged-mate regressions, negative-strand
coverage, resumability and input changes, default/sensitivity separation,
unavailable sources, failure states, ncORF reconstruction and PDF visual review.
