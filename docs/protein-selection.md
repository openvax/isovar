# Protein context selection

For each variant, Isovar can reconstruct several mutant protein candidates that
differ in length and read support. This page explains how much context it aims
for and which candidate it reports first. The policy does not change allele
counts, base-quality filters or assembly. All reported sequence comes from RNA:
missing context is never filled in from the reference.

## In short

- For peptides of length K (default 25), Isovar aims for **2K−1** residues
  (49 aa). That is enough to contain every K-mer overlapping a single changed
  residue.
- The default `balanced` preference considers candidates with at least **85%**
  of the best candidate's read support. It reports the one with the most
  mutation-containing K-mers.
- Every base of a reported cDNA must be covered by at least **2** reads.
- The result is shorter when coverage, a stop codon or the protein's end cuts it
  off. It can be shorter than K when no full K-mer meets the support budget, so
  check both length and mutation overlap before designing peptides.

## How much context

The target is 2K−1 residues: 15-mers → 29 aa, 25-mers → 49 aa, 30-mers → 59 aa.
With K−1 residues on each side of a single substituted residue, every
mutation-containing K-mer fits. `protein_sequence_length` sets an explicit
target instead.

Larger changes need more. A mutant interval W residues wide needs W + 2(K−1)
residues to contain every overlapping K-mer. A pure deletion has zero width, so
only windows that strictly span the deletion junction count: K−1 residues on
each side (48 for K = 25) give 24 junction-spanning 25-mers. The default target
is therefore a bound, not a promise of complete context for frameshifts or
multi-residue changes.

## Which candidate comes first

| Preference | Ranking |
|---|---|
| `balanced` (default) | Keep candidates with at least `min_protein_sequence_support_fraction` (0.85) of the best candidate's compatible read-name support. Among those, maximize mutation-containing K-mers. If none reaches K residues, prefer longer partial context. Ties go to support, then reads, mismatches and length. |
| `support` | Assemble one context length and rank by fragments, then reads, mismatches and length. With `protein_sequence_length=20` this is the policy from before Isovar 1.8. |
| `context` | Maximize mutation-containing K-mers with no support budget. It can pick weak or discordant RNA haplotypes, so it is not a shortcut to confidence. |

Candidates are ranked, never merged across conflicting haplotypes.
`max_protein_sequences_per_variant` (default 1; 0 keeps all) limits how many
are returned, after ranking. Allele read counts are reported separately from
each candidate's support.

## The two thresholds

1. **Support fraction (default 0.85).** A candidate's compatible read names,
   divided by the best candidate's, must reach this fraction.
2. **Coverage floor (default 2).** Every retained cDNA base must be covered by
   at least `min_variant_sequence_coverage` reads; 0 disables the floor. A
   variant that cannot meet it returns no mutant sequence, rather than one
   with the floor silently relaxed.

**85% is not 85% of per-base depth, or of all alt reads.** For example, 100
reads may support a short central region while only two extend into the flanks.
With a floor of five, those flanks are trimmed, even though all 100 names are
compatible with the longer candidate. Many short reads cannot rescue sparsely
covered flanks.

Neither measure counts independent molecules: paired mates may share a name, and
merging overlapping mates changes what is counted for coverage. When synonymous
RNA haplotypes give the same protein, their compatible names are combined, and
the coverage floor still applies to each RNA sequence. Compatible support can
include reads that span only part of a candidate, so it does not count
independently observed full-length peptides. Lowering the support fraction never
relaxes the coverage floor or the reference-matching filters.

## How candidates are generated

A longer context can lose support when noisy flanks fail to match the reference.
`balanced` and `context` therefore assemble a short ladder of context lengths: a
20-aa fallback, K, and at most six steps up to the target (20, 25, 29, 33, 37,
41, 45 and 49 aa for K = 25). Identical RNA sequences and identical
translations are merged, with their reads deduplicated. Each length requests
3 × length + 2 cDNA bases, enough for a partial leading codon and a centered
odd-length SNV context in every phase; `support` requests 3 × length + 3.

## Configuration

```python
from isovar import ProteinSequenceCreator, run_isovar

creator = ProteinSequenceCreator(
    protein_context_peptide_length=25,              # target: 49 aa
    protein_sequence_preference="balanced",         # or "support" / "context"
    min_protein_sequence_support_fraction=0.85,     # of compatible read names
    min_variant_sequence_coverage=2,                # reads at every cDNA base
)
results = run_isovar("variants.vcf", "rna.bam", protein_sequence_creator=creator)
```

The same options on the command line, with their defaults:

```text
--protein-context-peptide-length 25
--protein-sequence-preference balanced
--min-protein-sequence-support-fraction 0.85
--min-variant-sequence-coverage 2
```

For example:

- `--protein-context-peptide-length 30 --min-variant-sequence-coverage 5`
  targets 59 aa for 30-mers and requires five reads at every RNA base.
- `--min-protein-sequence-support-fraction 0.95` favors support more strongly.
- `--protein-sequence-length 20 --protein-sequence-preference support` restores
  the extraction and ranking from before Isovar 1.8.

Vaxrank passes its vaccine peptide length as `protein_context_peptide_length`
and lets Isovar derive the protein length. The exception is a Vaxrank
configuration with explicit padding, which asks for K + 2 × padding residues.

## Background

`balanced` has been the default since Isovar 1.8.0
([#219](https://github.com/openvax/isovar/pull/219)). The
[six-locus, two-source original-RNA audit](../tests/data/osteosarc/SAMPLE_AUDIT.md)
compares the three preferences against independent protein expectations.

This is a configurable sequence-selection policy, not an immunogenicity model.
The PGV-001 methods describe 25-mer design and RNA-based haplotype
reconstruction: [Rubinsteyn et al., 2018](https://www.frontiersin.org/journals/immunology/articles/10.3389/fimmu.2017.01807/full).
