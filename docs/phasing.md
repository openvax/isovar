# Phasing variants with RNA reads

Two variants are phased when their alt reads share at least
`min_shared_fragments_for_phasing` (default 2) fragments with compatible
alignments and pass the error-aware [cis test](#cis-or-trans) below. Both
partner lists and phase groups use all reference and alternate reads for that
test; protein-sequence links additionally need enough shared
protein-supporting fragments.

| `IsovarResult` property | Reads used |
|---|---|
| `phased_variants_in_supporting_reads` | All alt reads |
| `phased_variants_in_protein_sequence` | The reads behind the top protein sequence |
| `phase_group_from_supporting_reads`, `phase_group_from_protein_sequence` | As above, returning the connected `PhaseGroup`, which can include variants linked only through others |

A phase group is connected pairwise evidence, not one resolved haplotype.
`IsovarReadPhasing` and `IsovarMutantTranscript` expose these results through
Varcode's phasing and mutant-transcript interfaces.

## Which reads can link two variants

Complementary mates, variants on one spliced alignment, and supplementary
pieces whose reciprocal `SA` tags declare the same chimeric path can phase.
Matching read names in different read groups cannot.

## Cis or trans

Not being phased is not evidence of trans: Isovar only compares the variants in
its run.

`IsovarReadPhasing.in_cis(v1, v2)` answers from the fragments that cover both
loci. Each fragment carries one of four combinations: both alt alleles, either
alt allele alone, or neither. A fragment whose unmerged mates support
conflicting combinations is excluded from the phase counts, while its reads
remain in the result for inspection.

A combination counts when it is on at least `min_shared_fragments_for_phasing`
fragments, and on more fragments than reads showing the wrong allele would move
into it from its two neighbouring combinations. That second condition is a
one-sided binomial test at `phasing_error_rate` (1% by default; raise it for
ONT) and `max_p_value_for_phasing` (default 0.05).

`in_cis` returns:

- `True` when both alt alleles count together and not each alone;
- `False` (never on one molecule) when each counts alone and not together;
- `None` otherwise, including when all three count.

This is the four-gamete test (Hudson and Kaplan 1985) applied to read-phased
variant pairs, as Nik-Zainal et al. 2012 interpret them:

- A variant that arose in a subclone of the other's cells, on the other's
  copy, is cis with it on every molecule it's on. Cells without it still give
  the earlier variant's alt allele alone, so that alone isn't trans.
- Errors put a variant's alt allele on a small fraction of the other's alt
  fragments, which isn't cis.

## Error thresholds

`run_isovar`, `annotate_phased_variants` and `IsovarReadPhasing` accept
`phasing_error_rate` and `max_p_value_for_phasing`. The command-line options
are `--phasing-error-rate` and `--max-p-value-for-phasing`.

The run stores these thresholds in `result.phasing_settings`, and adapters
inherit them, including after cloning or pickling. Explicit adapter overrides
recompute partner lists. Legacy results without settings retain the default
thresholds and supplied lists. Results from runs with conflicting settings
require explicit overrides or reannotation.

### Per-locus calibration

For calibrated locus-specific errors, pass
`phasing_error_rates={variant: (ref_to_alt, alt_to_ref)}` to any of the three
entry points above. On the command line, pass
`--phasing-error-rates calibration.json`, a list of records:

```json
[{"contig": "1", "start": 123, "ref": "A", "alt": "G",
  "ref_to_alt": 0.10, "alt_to_ref": 0.02}]
```

These numbers illustrate the format; they are not recommended platform rates.
Coordinates and alleles must match the normalized input variants (1-based
start).

Both directions use at least the scalar `phasing_error_rate` floor. Each allele
combination is tested against the larger of its two incoming error rates; this
conservatively bounds leakage under the independent-error model. An alt-to-ref
bias therefore does not inflate the error threshold for both-alt support. This
calibration currently applies to pairs of variants in the run; matched-germline
majority decisions keep their separate interpretation
([below](#matched-germline-variants)).

Retain the control dataset, library, basecaller, context and calibration method
alongside the profile. For mixed read groups, supply a bound covering every
included group or analyze them separately.

### Limits of the error model

The default 1% is not an automatic ONT or PacBio calibration. Non-ref/non-alt
reads alone cannot identify target-specific A→G editing or a homopolymer
deletion bias. Nanopore error rates are directional and depend on sequence
context ([primary study](https://doi.org/10.1186/s12864-024-10440-w)).

Use appropriate control-derived bounds; do not interpret an uncalibrated
definite call as evidence that these errors were excluded. Correlated RT/PCR
artifacts and uncertain calibration are not covered by the independent
binomial model.

## Matched germline variants

With `run_isovar(germline_variants=...)` (or `--germline-vcf`), Isovar also
collects reads at each matched germline variant that a variant's alt reads
cover, outside skipped introns. They are in `IsovarResult.germline_read_evidence`.

Only fragments carrying the somatic alt allele count: with the germline alt
allele they are cis, with its reference allele trans. Fragments with the
somatic reference allele say nothing about its copy, because the germline alt
allele also comes from cells without the somatic variant (normal
contamination, other subclones) and from both copies of a homozygous variant.
A call needs the cis or trans count to reach
`min_shared_fragments_for_phasing` and exceed the other.

Without decisive reads, a germline edit in the variant's top assembled protein
is cis. If the reads say trans but the protein has the edit, the answer is
unknown (`None`).

Unmerged mates also nominate germline loci through their primary aligned blocks
(MC or a matching mate lookup), excluding skipped introns and unsequenced
inserts. Varcode's `MolecularPhaseResolver` uses this method for cis/trans
decisions.

The same matched normal also labels edits in the assembled RNA as germline;
see [other variants in the assembled RNA](assembly-and-translation.md#other-variants-in-the-assembled-rna).
