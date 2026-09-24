# Protein context selection

This policy chooses which RNA-derived protein context to report for each
variant. It does not change allele counts, base-quality filters or assembly
defaults. All sequence is RNA-derived; the reference never pads missing context.
It has been the default since Isovar 1.8.0 ([#219](https://github.com/openvax/isovar/pull/219)).

## Context target

For peptide size K (default 25), the target length is **2*K-1**: 15mers target
29 aa, 25mers 49 aa and 30mers 59 aa. For a single substituted residue, K-1
residues on each side include every mutation-containing K-mer. An explicit
`protein_sequence_length` overrides the derived target.

Only windows that intersect the mutant interval count, or, for a zero-width
deletion junction, windows that strictly span it. A mutant interval of width W
needs W + 2*(K-1) residues for every overlapping K-mer, and a pure deletion
junction needs K-1 residues on each side (48 for K=25), giving 24
junction-spanning 25mers. The default target is therefore bounded, not a
promise of complete context for every frameshift or multi-residue change. Stops,
protein boundaries and RNA coverage can all give shorter output.

## Preferences

- **`balanced`** (default): among mutant candidates with at least
  `min_protein_sequence_support_fraction` (default 0.85) of the best candidate's
  compatible read-name support, maximize the number of mutation-containing peptide
  windows. If none reaches K residues, prefer longer partial context within the
  same budget. Support, reads, mismatches and length break ties.
- **`support`**: assemble one context length and rank by fragments, reads,
  mismatches and length. With `protein_sequence_length=20` this reproduces the
  pre-1.8 policy.
- **`context`**: maximize peptide windows without the support budget. It can
  select weak or discordant RNA haplotypes and is not an evidence-confidence shortcut.

Independently, every retained cDNA base must be covered by at least
`min_variant_sequence_coverage` reads (default 2; 0 disables the floor). A large
number of short compatible reads cannot rescue sparsely covered flanks, and a
variant that cannot meet the floor returns no mutant sequence rather than a
silently relaxed one.

Candidates are ranked, never merged across conflicting haplotypes, and the
`max_protein_sequences_per_variant` cap (default 1; 0 keeps all) applies last.
Allele read counts are reported separately from candidate support.

## Candidate generation

Longer context alone can lose reference-matching support when noisy flanks fail
to match. The `balanced` and `context` preferences therefore assemble a bounded
ladder of context lengths: a 20-aa fallback, the peptide length, and at most six
increments up to the target (20, 25, 29, 33, 37, 41, 45, 49 for K=25). Identical
RNA sequences and identical translations are merged, and their original reads
deduplicated. Each length requests 3 × length + 2 cDNA bases, enough for a partial
leading codon and a centered odd-length SNV context in every phase; `support`
requests 3 × length + 3.

## Configuration

```python
from isovar import ProteinSequenceCreator, run_isovar

creator = ProteinSequenceCreator(
    protein_context_peptide_length=25,              # target: 49 aa
    protein_sequence_preference="balanced",         # or "support" / "context"
    min_protein_sequence_support_fraction=0.85,     # of compatible read names
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

For example, `--protein-context-peptide-length 30 --min-variant-sequence-coverage 5`
targets 59 aa for 30mers while requiring five reads at every retained RNA base.
Use a fraction of 0.90 or 0.95 to favor support more strongly.
`--protein-sequence-length 20 --protein-sequence-preference support` restores
the pre-1.8 extraction and ranking.

The top result may still be shorter than K when no complete window meets the
support budget, so downstream design must check both length and mutation overlap.
Lowering the support fraction never relaxes the coverage floor or the
reference-matching filters. Compatible support can include reads spanning only
part of a candidate, so it is not a count of independently observed full-length peptides.

Vaxrank passes its vaccine peptide length as `protein_context_peptide_length`
and leaves the protein length to this derivation, unless its configuration
sets an explicit padding (K + 2 × padding residues).

## What the two thresholds mean

1. Candidate-compatible read names divided by the best candidate's compatible
   read names must be at least the support fraction.
2. Every retained cDNA base must meet the absolute coverage floor.

**85% does not mean 85% of per-base depth or of all alternate reads.** For
example, 100 reads can support a short central region while only two extend into
the flanks. A floor of five trims those flanks even though all 100 names are
compatible with the longer candidate. Paired mates may share a name, and
overlapping-mate merging changes the objects counted for coverage; neither
measure is an independent molecule count. When synonymous coding-RNA haplotypes
group into one protein, compatible names are unioned, and the coverage floor
still applies to each retained RNA sequence.

The [six-locus, two-source original-RNA audit](../tests/data/osteosarc/SAMPLE_AUDIT.md)
compares these preferences against independent protein expectations.

This is a configurable sequence-selection policy, not an immunogenicity model.
The PGV-001 methods describe 25-mer design and RNA-based haplotype reconstruction:
[Rubinsteyn et al., 2018](https://www.frontiersin.org/journals/immunology/articles/10.3389/fimmu.2017.01807/full).
