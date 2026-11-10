# One effective code for native ORF reconstruction (#700)

Use the existing `ReferenceContext.translation_table` for candidate, focal-allele
counterfactual, and reference-window decoding. Do not add another biological
code resolver. Preserve `reference.translation_table` as source metadata,
including undeclared None. The sequence-only helper keeps its table-1 default.

Emit `isovar.variant_rna_orfs.v2`. Each retained protein candidate records the
effective integer `translation_table` and `translation_table_source`:
`orf_model` for an explicit source declaration, otherwise `reference_context`.
Include the effective code in candidate identity. Unresolved models do not gain
invented protein candidates or asserted source codes.

Expose a small public `variant_orf_translation_table(candidate, reference,
schema=...)` reader. V2 requires the code and consistent provenance, rejects
unsupported codes and conflicts with explicit model declarations. V1 replay
retains its original nullable source-code convention; it is not reinterpreted.
Use this reader in both Isovar's ORFdata snapshot adapter and Vaxrank's native
adapter, including frame grouping, validation and stored candidate metadata.
Vaxrank decoding should call the existing shared translation helper directly.

Validate native MT-ND5 (+), MT-ND6 (-), and TP53 against sequence-validated
explicit-code controls; candidate proteins, reference and counterfactual
windows, and RNA support must agree. Test explicit source codes 2 and 3,
nuclear NuMT-like carriers, invalid/missing export codes, provenance conflicts,
v1 replay, JSON round trips, ORFdata snapshots and Vaxrank validation/replay.

A nuclear insertion of mitochondrial-like sequence does not become a
mitochondrial transcript: mammalian nuclear carriers use table 1 unless an
explicit validated source model says otherwise. This is a compartment/code
regression, not a NuMT detector. Existing osteosarc diagnostics report that
NuMT origin is not ruled out; general competitive alignment is outside this
fix. Do not claim mitochondrial origin from MAPQ, NH, gene name or similarity.

Release as Isovar 1.89.0 with a coordinated Vaxrank release requiring it.
Both branches must pass lint, full tests and CI, then merge and deploy.

Primary sources: [NCBI codes](https://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi),
[Maude et al. 2019](https://doi.org/10.3389/fcell.2019.00201),
[Laricchia et al. 2022](https://pmc.ncbi.nlm.nih.gov/articles/PMC8896463/).
