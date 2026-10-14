# Partial-read support scoring (#26)

## Scope

Add an opt-in, catalog-relative fractional fragment support score to protein
hypothesis exports and pooled/cell attribution. Retain raw counts, reconstruction,
existing ranking and compatibility assignments. The score is an evidence summary,
not an abundance estimate, probability, molecule count or adaptive pooling policy.

Use the existing anchored edit/transcript-path comparison for bulk and cell
outputs. Move that comparison into a shared module without changing its default
behavior. All candidates within the explicit edit tolerance remain alternatives;
a lower edit count alone does not manufacture certainty about a noisy base.

## Score and counting unit

- Re-evaluate original ref/alt/other observations against the whole discovered
  candidate catalog, not just reads used to assemble each candidate.
- Group by source-scoped (read group, QNAME) fragment. Mates, merged mates,
  repeated records and alternative placements must never create extra units.
- For each observation, union candidates across alternative placements. If an
  alternative placement disagrees on the compatible candidate set, retain that uncertainty
  and withhold the fragment's score. Reads with no usable comparison are reported
  separately, not fabricated as support. Joint mate support is the intersection
  of their candidate sets; conflicting mates cannot create a chimera.
- Each scored fragment contributes 1/K to each of K compatible RNA candidates.
  For proteins, first project the joint RNA candidate set to distinct protein
  hypothesis IDs, then divide by that count. Synonymous transcripts cannot
  inflate a protein's support. Report unique and ambiguous supporting fragment
  counts alongside the fractional sum and all underlying fragment assignments.
- Unique means unique within the exported catalog and chosen compatibility
  settings. A tiny catalog or truncated candidate cap cannot establish global
  specificity; export catalog completeness.
- Do not divide by full candidate length. Unobserved extensions are neutral,
  and a short discriminating read can contribute one unit. Report compared
  intervals and edit counts already produced by the shared comparator.
- Optional focal-allele quality annotations stratify the same assigned score:
  any failed source read fails the fragment category; otherwise any unassessed
  source read makes it unassessed. Do not score flanking differences with focal
  quality, substitute qualities or introduce platform defaults.

## Integration

Add `partial_read_support=True` to protein and cell export APIs and
`--partial-read-support` to `protein-hypotheses`. Shared comparison settings
have CLI aliases `--support-max-edits` and `--support-min-overlap` (existing
cell options remain valid). JSON contains the scored catalog, original fragment
assignments, summaries and a separate dense score rank; raw Isovar ordering stays
unchanged. TSV gains score/rank columns only when opted in. Cell exports retain
sample/header boundaries and summarize scored fragments within each resolved
cell. Cells with conflicting/missing labels remain explicitly unresolved.

## Validation

Test short informative versus shared reads; score conservation; observed-only
windows; missing QUAL/platform invariance; read order, duplicates, mates and
alternative placements; reference competitors; conflicting or insufficient
observations; synonymous RNA/protein projection; incomplete/empty catalogs;
both strands; cell scopes; quality categories; CLI/TSV opt-in and default parity.
Use the existing offline public RNA and pinned transcript fixtures to exercise
real reads through the shared scorer.

## Scientific context

[RSEM](https://pmc.ncbi.nlm.nih.gov/articles/PMC3163565/) models ambiguous
read assignments, fragment generation and abundance jointly. This PR deliberately
implements a transparent equal allocation over compatible hypotheses, without
claiming the posterior probabilities or abundance estimates of that model.
[The kallisto paper](https://doi.org/10.1038/nbt.3519) likewise separates read
compatibility from the subsequent probabilistic quantification step. Calibration,
error-weighted likelihoods and adaptive pooling remain separate future work.
