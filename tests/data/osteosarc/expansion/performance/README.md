# High-depth RNA performance (#227)

## Spec

Profile the checksum-verified full regional osteosarc sources, not the bounded
offline fixtures. Measure production reconstruction and independent protein
validation separately, including wall time, peak process RSS, and Python
function profiles. Keep timeouts explicit and retain completed evidence.

Optimize measured hotspots without changing the candidate set, original read
objects/multiplicity, coverage, translations, or protein ordering. Compatibility
must remain anchored at the allele and non-transitive: a short shared core does
not connect conflicting longer haplotypes. No downsampling, changed thresholds,
reference padding, or new biological assumptions.

Verify against the exhaustive matcher using randomized and adversarial cases,
and compare complete outputs on original reads. Include empty flanks, deletions,
duplicate names, conflicting branches, coverage trimming, and adaptive lengths.
Benchmark full-depth mitochondrial sources and smaller nuclear controls. Rerun
both modes of the five interrupted #218 mitochondrial cells, with independently
checked proteins and provenance-bearing results. Document remaining limits.

Run lint, the complete offline test suite, and GitHub CI; include a patch version
bump and link the PR to #227 and #218.

## Profile-driven refinement

The first indexed run moved the bottleneck into per-read coverage updates.
Accumulate exact coverage from the disjoint original-read groups while matching,
and carry an independent copy of coverage through interval trimming. Preserve
the per-read interval calculation for standalone `VariantSequence` objects.
Disabled logging must not format coverage arrays eagerly.

The remaining collection-stage memory peak is tracked separately in #229.
Some exploratory runs overlapped other workstation jobs; retain timed-out
measurements and report both process CPU time (where recorded) and wall time.
