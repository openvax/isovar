# Reproducible CTA target selection

Issue #641: explicit genes are currently required and annotation preparation is
embedded in the CLI. Introduce a public target/reference preparation function
used by API, CLI, and cohort manifests. The named `oncoref-cta` selection reads
the published OncoRef species-scoped CTA decisions, retaining every nominated
gene locus and every available transcript in the supplied annotation.

Record provider version, policy, snapshot, source inventory, gene mapping,
unavailable models, and reference checksums. Never reinterpret a missing stable
gene ID by a similar symbol. Explicit user-supplied symbols remain supported.
Selection is research nomination, not Vaxrank admission or translation proof.

Keep explicit target-model input compatible. Core dependency adoption must be
coordinated with Vaxrank/PirlyGenes exact OncoRef pins before release. Validate
API/CLI/cohort parity, missing annotation coverage, reference identity, and
receipt invalidation. Per-gene discovery budgets and scalable matching are
required before the full catalogue is used for the sample survey.
