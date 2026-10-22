# Multisource evidence plots

## Denominators

Keep the existing VCF small-variant funnel. Add an exon-involved subset of the
ALT-read histogram, with and without zero-support records. Its denominator is
all small-variant source alleles whose changed reference span intersects any
Ensembl 87 exon; insertion flanks remain candidate overlaps. Do not label a
symbolic SV as an ALT-spanning small-variant read.

For nominated SV source records, annotate any involved exon using the interval
of an intrachromosomal DEL/DUP/INV, or the nominated breakpoint positions for
INS/TRA/BND. Plot split-alignment candidate counts by SV class, both including
zero and conditional on a positive count. This is a screening metric, not a
validated RNA junction or productive SV protein. Keep SV source sets distinct.

For viral exports, plot original viral-contig aligned and sequence-assigned
fragment counts by contig in their own panel. An absent viral export is
"not assessed"; a zero within a supplied export is an assessed zero for the
selected contig/reference only. Do not apply a human-exon filter or an ALT-read
metric to viral paths.

## Report composition

The existing `isovar evidence-report` retains its VCF/BAM workflow and adds the
exon-involved plot. A new `isovar antigen-report` composes its per-call TSV with
optional versioned Isovar SV, fusion, viral and antigen-isoform exports. It
validates schema and sample identity, records source files, and shows each
category separately. Counts for nominated DNA events, supplied fusion RNAs,
selected viral contigs and targeted genes are never summed into a single
"total antigens" number. Reference-only isoforms are not counted as observed
tumor products. ncORF candidate overlaps and exploratory SV ORFs remain
separate from demonstrated protein expression.

## Verification

Test interval/breakpoint exon semantics, plot denominators, schema/sample
validation, viral fragment units, and one synthetic mixed-source PDF. Rebuild
both HCC1395/COLO829 PDFs, inspect every page, then run lint, tests and CI.
