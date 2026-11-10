# Standard genetic-code default

Restore `translate_cdna(sequence)` and its existing positional start flag by
defaulting the keyword-only `translation_table` to 1. This is the standard
code used for mammalian nuclear translation. Explicit None retains its existing
table-1 behavior; explicit alternative codes still use the same Varcode NCBI
registry. Do not restore duplicate tables or the retired mitochondrial flag.

Test omission against independent NCBI assignments for all 64 codons, standard
initiation versus elongation, and an explicit table-2 override. Keep invalid-code
validation. Update the public guide and release as a patch correction to 1.88.0.

Separately verify #700 with native Ensembl annotated nuclear and mitochondrial
models and the public RNA reconstruction API, on both 1.87.1 and 1.88.0. The
sequence-only helper's standard default cannot infer transcript compartment;
native reconstruction must use its available reference context. Preserve
undeclared source metadata rather than treating this API-default correction as
a fix for effective-code resolution. Record any factual corrections on #700.

Primary source: [NCBI genetic-code definitions](https://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi).
