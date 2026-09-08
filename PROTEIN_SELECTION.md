# RNA-supported vaccine context

Proposed in Isovar 1.8.0 / PR #219 (not yet released). This changes protein
context generation and selection, not allele counts, base-quality filters,
RNA assembly defaults or clinical scoring. All sequence remains RNA-derived;
there is no reference padding. The PR is left open for review.

## Default and alternatives

- The context target is **derived from desired peptide size K: 2*K-1**, not
  fixed at 49. Thus 15mers target 29 aa, 25mers target 49 aa, and 30mers
  target 59 aa. K defaults to 25; actual output adapts to RNA support and
  per-base coverage, and may be shorter than either target or peptide size.
  For a single substituted residue, 24 aa on each side includes all 25
  mutation-containing 25-mers. Length without mutation position is not
  sufficient. Count only windows intersecting the mutant interval, or
  strictly spanning a zero-width deletion junction. Stop/protein boundaries
  can make fewer windows genuinely available. Wider mutations may require
  a larger explicitly requested context to cover every placement: a mutant
  interval of width W needs W + 2*(K-1) residues for all overlapping Kmers.
  A pure deletion junction needs K-1 residues on each side (48 for K=25),
  giving 24 junction-spanning 25mers, not 25. The default 49 is a bounded
  target, not a promise of complete context for every frameshift or MNV.
- Default preference: `balanced`. Among mutant candidates retaining at
  least 85% of the best candidate's compatible read-name support, maximize
  the number of mutation-containing peptide windows. If none reaches the
  requested peptide length, favor useful partial context within that same
  support budget; never claim a full vaccine window or relax the budget.
  Use existing support/mismatch ordering to break ties. The fraction is a
  configurable engineering tolerance, not a biological confidence estimate.
- Independently require **at least 2 RNA read objects at every retained
  cDNA base** by default, using `min_variant_sequence_coverage`. This is a
  configurable absolute coverage floor, not a fraction of candidate support.
  A high count of short compatible reads cannot rescue sparsely covered
  flanks below this floor. If the variant cannot meet the floor, return no
  usable mutant sequence rather than silently lowering it.
- Alternatives: `support` preserves the existing support-first single-scale
  algorithm (explicit length 20 reproduces the prior policy); `context`
  prioritizes useful peptide windows without the relative support budget.
  Keep all candidates internally and apply the existing output cap last.
- Preserve allele read counts separately from candidate support. Never add
  incompatible flanks or claim every compatible name spans a whole peptide.
  No abundance/probability reweighting is introduced by this selection change.

## Candidate generation

Increasing context alone can lose reference-matching RNA support. For the
new preferences retain a bounded ladder of shorter RNA contexts along with
the target: a 20-aa fallback, the peptide length, and at most six increments
up to the target. For 25/49 this is 20,25,29,33,37,41,45,49. Merge only
identical RNA sequences and identical translations, deduplicating original
read objects/names. Do not assemble across conflicting haplotypes.

Use two extra nucleotides to cover partial leading codons for new modes;
this allows a fully covered odd-length SNV context to center the changed
codon in every phase. Keep the old extraction budget in `support` mode.

## Configuration

```python
from isovar import ProteinSequenceCreator, run_isovar

creator = ProteinSequenceCreator(
    protein_context_peptide_length=25,              # auto target: 49 aa
    protein_sequence_preference="balanced",         # or support / context
    min_protein_sequence_support_fraction=0.85,      # compatible names
    min_variant_sequence_coverage=2,                # absolute per-base floor
)
results = run_isovar("variants.vcf", "rna.bam", protein_sequence_creator=creator)
```

Equivalent CLI options (defaults shown):

```text
--protein-context-peptide-length 25
--protein-sequence-preference balanced
--min-protein-sequence-support-fraction 0.85
--min-variant-sequence-coverage 2
```

For example, use `--protein-context-peptide-length 30
--min-protein-sequence-support-fraction 0.85 --min-variant-sequence-coverage 5`
to target 59 aa for 30mers, retain 85% of candidate support, and stay above
5-read coverage at every retained RNA base. Use fraction 0.90 or 0.95 to
prioritize relative retention more strongly. `--protein-sequence-preference context` removes the relative
budget; it can select weak, discordant RNA haplotypes and is not recommended
as an evidence-confidence shortcut. `--protein-sequence-length 20
--protein-sequence-preference support` restores the historical extraction
and support-first ranking. For another peptide size K, omission of an
explicit protein length derives 2*K-1 automatically. An explicit
`--protein-sequence-length` always takes precedence. Increase
`--max-protein-sequences-per-variant` to inspect alternatives; 0 returns all.

The default top-one result may still be shorter than K if no complete
window meets the support budget. Downstream design must check both length
and mutation overlap. Lowering the relative budget is explicit; it never
relaxes minimum coverage or reference-matching filters. Compatible support
can include reads spanning only part of a candidate, so it is not a count
of independently observed full-length vaccine peptides.

### What the two thresholds mean

These are two separate requirements:

1. Candidate-compatible read names / best candidate-compatible read names
   must be at least the configured fraction (default 0.85).
2. Every retained cDNA base must meet the absolute read-object coverage floor
   (default 2; explicitly set 0 to disable that floor).

**85% does not mean 85% of per-base depth or of all alternate reads.**
For example, 100 reads can support a short central region while only two
extend into the flanks. A floor of five trims those flanks even though all
100 names are compatible with the longer candidate. Paired mates may share
a name, and overlapping-mate merging changes the objects counted for
coverage. Neither measure is automatically an independent molecule count.
When synonymous coding-RNA haplotypes group into one protein, compatible
names are unioned; the per-base floor still applies to each retained RNA
sequence. The audit records both measures explicitly.

## Integration and verification

- Validated Python/CLI options control preference, peptide size, and minimum
  relative support. Explicit protein length is respected. Use vaxrank's configured
  vaccine peptide size when its argparse namespace supplies one and no
  Isovar-specific peptide-size override is given.
- Vaxrank currently explicitly requests vaccine length + 2*padding (35 aa
  for its defaults). Do not silently override that explicit caller value;
  see [vaxrank #415](https://github.com/openvax/vaxrank/issues/415). Padding 12
  already requests 49 aa for 25-mers without a vaxrank code change.
- Minor version 1.8.0 signals the intentional default behavior change.
- Exhaustive window-count oracle tests, thresholds/invalid inputs/order
  invariance, zero-width deletions, frameshifts, stops, all three SNV phases,
  CLI propagation and historical support-policy regression coverage.
- The [six-locus/two-source original-RNA audit](tests/data/osteosarc/SAMPLE_AUDIT.md)
  uses independent protein expectations and records support/context tradeoffs
  and no-alt/error states. The previous 20-aa comparison is explicitly named
  historical policy, not current defaults. Full-source checksums and an
  independent regeneration verify the report outside offline CI.

Scientific scope: this is a configurable sequence-selection policy, not an
optimized immunogenicity model. The PGV-001 methods describe 25-mer design
and RNA-based haplotype reconstruction:
https://www.frontiersin.org/journals/immunology/articles/10.3389/fimmu.2017.01807/full
