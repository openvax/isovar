# Separate RNA evidence figures

The two cell-line PDFs will show separate nested funnels for small-variant
alleles and nominated structural-variant source records. Each percentage uses
its own DNA input count. The small-variant path ends with an observed RNA ALT
read and a local protein context. The SV path distinguishes RNA coverage
near a nominated breakend, a same-read split or cross-end mate candidate, a
resolved RNA junction, and a translated SV product. Coverage alone is not
event support; a candidate does not establish the junction sequence.

Read-support histograms will put each variant class and its count in a fixed
column with matching colors, eliminating position-dependent slash labels.
The first column remains the read-support bin. Color is redundant with the
printed class name and count.

Five high RNA-support small-variant examples will show a pinned Ensembl 87
reference protein context and the top Isovar mutant local context, aligned at
the mutation boundary. Shared flanks and changed amino acids will have
different colors. For frameshifts the downstream columns are a boundary-
anchored comparison, not a homology alignment or full-length translation.
Each example must name the transcript and retain the original-read support
count. A fixed reference-context table keeps PDF regeneration offline.

Verification: synthetic funnel and figure tests, source-table count assertions,
rendered-PDF text checks, visual inspection of both PDFs, lint and full tests.
