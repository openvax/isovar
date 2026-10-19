# Isovar #522: retain and compare alternative consequences

## Contract

- Keep Varcode outcome-set objects and their provenance through coding and
  sequence-availability filters. Coding eligibility uses the supported
  `modifies_coding_sequence` contract, including potentially coding splice sets
  and unresolved status; nested alternatives must not be hidden by a primary.
- A predicted protein is available when any candidate has a sequence. Traverse
  `MultiOutcomeEffect.candidates`, including concrete self-candidates, safely.
- Compare the RNA-derived protein against each candidate using its own protein
  and mutation offsets. A matching candidate establishes compatibility, not
  proof of its mechanism. No match with an unresolved candidate remains unknown.
- Existing scalar trimmed fields return a value only when every candidate has
  the same available value. Preserve ordinary single-effect behavior. Keep the
  full original effect on the result; do not replace it with a selected candidate.

## Verification and release

Reproduce the SCP2 example from #522 with Ensembl 75; cover coding/noncoding,
silent and unresolved candidates, nested/self candidate sets, multiple known
proteins, candidate-specific offsets, exact/no/unknown matches and ordinary
effects. Run focused tests, downstream Vaxrank integration, lint and the full
suite. Bump the patch version with release notes, obtain green CI, merge, deploy
from clean master, verify PyPI, then update the shared installation normally.

## Sources

Varcode's `MultiOutcomeEffect`, `EffectCandidate`, `SpliceOutcomeSet` and
`SpliceMechanismEffect` public contracts define the implementation.
The reason to preserve alternatives and uncertainty is consistent with Dawes
et al., *SpliceVault predicts the precise nature of variant-associated
mis-splicing*, Nature Genetics 55, 324–332 (2023),
https://doi.org/10.1038/s41588-022-01293-8 (multiple observed outcomes per variant).
No new splice prediction model is introduced.

## Remaining authorized maintenance

After Isovar, resolve PyEnsembl #470 against the in-flight cache compatibility
PR #474 without changing existing checkouts. Then fix Varcode #604 with
workflow/PR-scoped cancellation, unique groups for main pushes, and a live
two-commit cancellation check. Each repository ships independently.
