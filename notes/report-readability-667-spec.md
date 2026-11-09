# Reporting on shared producer contracts (#667)

Implement after #681 ships. Preserve the native RNA measurement subject and source; do not invent one interchangeable antigen count.

- Add named call/gene/effect rows to evidence PDF/TSV and antigen reports; preserve measured zeros versus unavailable values. Avoid duplicate all/exonic pages when the populations are identical.
- Make antigen reports source-aware: repeat small-variant inputs, explicitly map an input export to its input_sample_id while keeping the report sample label separate, retaining each source/library's counts. Reuse cohort source identity rules; do not create another manifest/sidecar format. Keep single-sample CLI behavior compatible.
- Include per-event/path/junction rows for SV/fusion and named CTA/viral paths with original reads/fragments, compatibility/uniqueness and full-span evidence where present. Use producer fields and public support contracts; neither pool overlapping paths nor count assembly voters as complete-input support. Retain all rows in JSON and paginate full tables without arbitrary twelve-row truncation.
- Use flowing PDF sections instead of forced per-category pages. Render representative sparse and dense reports, inspect every page, and test stable bytes/relocatable artifact listings.
- Plot transcripts by RNA support then focal exon overlap; show REF/ALT/OTHER denominators and explicit producer no-protein reasons. Label protein groups using covered_ranks. Add an explicit literal-variant selector for multi-variant inputs, preserving existing single-variant behavior.
- Improve count-axis tick labels and bar values, fusion depth traces and external legends; condition intronic explanatory text on the model's actual introns. Set ReportLab invariant per document and emit report_files.json paths relative to its root.

Some cohort axis handling is already present (_set_count_ticks for positive-depth survival curves); improve only remaining unreadable panels. Scientific predicates/normalization do not change. Primary-source review is required for any discovered scientific-semantic change, which should be filed separately rather than folded into reporting.

Verification: regressions for source mapping/rejection, named rows/count units/zeros, multi-variant selection, row ordering and group ranks, deterministic rendering and relative paths. Existing producer, CLI I/O and report tests plus lint/fullsuite/CI. Merge and PyPI release required by AGENTS.
