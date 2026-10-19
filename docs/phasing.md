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

#### Estimate a profile from control BAMs

`isovar calibrate-phasing --manifest controls.json --output profile.json
--report calibration-report.json` collects the same allele evidence as Isovar
from independently truth-labeled, indexed control BAMs. It counts each physical
fragment once, excludes fragments with conflicting or other-allele calls from
the binomial denominator, and reports those exclusions. It computes separate
one-sided exact [Clopper–Pearson](https://doi.org/10.1093/biomet/26.4.404)
upper bounds for ref→alt and alt→ref errors. The confidence level defaults to
0.95; a Bonferroni correction covers both directions, every included read
group and every control source at each locus. The emitted locus profile takes
the largest bound across sources and groups and works with `--phasing-error-rates profile.json` or
`calibrate_phasing_error_rates(...).rates` in Python.

The manifest names each target BAM and every control source. For example:

```json
{
  "confidence": 0.95,
  "reference_fasta": "GRCh38.fa",
  "targets": [{
    "variant": {"contig": "1", "start": 123, "ref": "A", "alt": "G"},
    "context_start": 120,
    "target_bam": "target.bam",
    "read_groups": {
      "rg1": {"platform": "ONT", "basecaller": "dorado-1",
              "library": "RNA-protocol-1", "biological_context": "matched-tissue",
              "reference": "GRCh38", "context": "TACAGTA"}
    },
    "controls": [
      {"truth": "ref", "source": "known-ref-sample", "bam": "ref.bam",
       "read_group": "rg1", "platform": "ONT", "basecaller": "dorado-1",
       "library": "RNA-protocol-1", "biological_context": "matched-tissue",
       "reference": "GRCh38", "context": "TACAGTA"},
      {"truth": "alt", "source": "known-alt-sample", "bam": "alt.bam",
       "read_group": "rg1", "platform": "ONT", "basecaller": "dorado-1",
       "library": "RNA-protocol-1", "biological_context": "matched-tissue",
       "reference": "GRCh38", "context": "TACAGTA"}
    ]
  }]
}
```

Paths are relative to the manifest. The FASTA must be indexed, and each
`context_start` is the 1-based position of the first base in the declared
sequence context. The command checks the context and reference allele against
the FASTA. Every target read group needs both true-ref and true-alt controls.
Read-group IDs and platform tags must match BAM headers; contig lengths and
available MD5 digests must agree. The manifest's basecaller, library,
biological context and reference label are declarations: verify them against
your experiment before use. An A→G editing
control needs matched biological conditions; a homopolymer indel needs the
same local sequence context and library/aligner behavior. The command rejects
missing truth, missing groups and missing callable coverage; conflicting
fragment observations are reported and excluded. No rate is inferred from tumor allele fractions or
non-ref/non-alt reads. The report records counts, exclusions, source BAM paths,
the manifest digest, metadata, per-group bounds and confidence. A control with 100% wrong callable
fragments cannot supply a rate below one and is rejected.

The command uses Isovar's default read-collection filters. If the target run
changes duplicate, secondary-alignment, mapping-quality, soft-clip, mate-merge
or missing-base-quality filters, set the matching `read_collector_options`
object at the top level of the manifest. The report records the supplied
options; use the same options when running Isovar. Other custom read filters
need calibration through the Python API with evidence collected using that
filter.

The binomial bound assumes independent fragments and controls representative
of the target. PCR/RT chimeras, repeated molecules, sample-specific editing
and overdispersion can violate that assumption. A control-derived profile is
evidence for the stated controls and conditions, not a guarantee for other
samples. Primary studies describe [directional, context-dependent nanopore
errors](https://doi.org/10.1186/s12864-024-10440-w) and
[tissue-dependent RNA editing](https://rnajournal.cshlp.org/content/19/2/257.full.html).

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
