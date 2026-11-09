# Shared reconstruction and support contracts (#681, then #667)

## Outcome

Small literal mutations, structural adjacencies, CTA isoforms and viral paths
should use public shared RNA machinery. Event-specific adapters select seeds,
references and interpretations; they must not maintain separate copies of read
identity, placement geometry, exact path assembly or support accounting.
Reports consume producer measurements rather than reconstructing counts.

## Delivery order

1. Ship #681: extract the shared engine from `sv_rna`, publish documented
   discovery/support contracts, migrate all package consumers, and eliminate
   cross-module imports of private package names. Keep existing public entry
   points and serialized evidence semantics stable.
2. Ship #667 on that foundation: consistent named antigen rows and explicitly
   labelled support units in reports, plus the issue's plot/CLI usability fixes.
   Detailed CTA/viral paths must remain inspectable instead of aggregate-only
   summaries. Reuse the existing original-witness support records.

PR #693 already owns #692's common local-feature count and experimental
normalization API. Do not build a competing count schema or normalization model.
Inspect its integration as it lands, and migrate any shared-engine imports with
the other consumers. #689 and #690 are separate correctness/performance reports;
retain their cases and boundaries without disguising behavioral changes as a
mechanical rename.

## Engine boundaries

- Original-library collection/index ownership and read policy remain shared.
- Move observation geometry, placement/run indexing and exact overlap extension
  out of the SV-specific module into cohesive public engine modules.
- Keep SV adjacency classification, seeding and result interpretation in SV
  adapters. Reuse the existing transcript-path frame layer.
- Public discovery and complete-input support functions serve literal-allele,
  CTA and viral adapters. The ordinary allele-anchored small-variant route keeps
  its explicit locus/coding constraints; use the same identity, support and
  geometry contracts where applicable, without inventing an SV for an SNV.
- Give already-shared auxiliary helpers meaningful public names in the module
  owning their concept. Rename shared audit/file-identity modules. Do not add
  compatibility aliases for previously private cross-module names.
- Add an architectural check rejecting new cross-module private imports, with
  independent public-consumer regressions for the engine contracts.

## Scientific invariants

Preserve sample/library-scoped original identities; alternative placements of
one segment cannot supply independent assembly witnesses. Keep unique,
compatible and unassigned fragments separate. Discovery budgets do not truncate
complete-input support. Unmeasured and unavailable results are not zeros.
Contributors, per-base depth, allele/junction support and complete-sequence
witnesses remain distinct measurements. Paired unsequenced gaps do not donate
sequence or become single-read spans. Shared alternative paths are not additive.

This first PR changes ownership and contracts, not scientific predicates or
normalization. Any discovered semantic defect needs an explicit issue and
regression, and primary-source review before a scientific change.

## Verification and release

Record the existing consumer/dependency map before edits. Exercise public
observation, placement, assembly and support APIs with independent fixtures,
including both strands, ambiguous placements, overlap assembly and incomplete
support. Run relevant existing literal-allele, SV bounded/dense, CTA, viral,
audit and small-variant regressions; run `./lint.sh` and `./test.sh` before merge.
Check downstream public imports and serialized schema compatibility. Reserve
the release version against current master/PRs before publishing. Each PR must
pass GitHub CI, merge, and deploy from clean master, with PyPI/tag verification.
