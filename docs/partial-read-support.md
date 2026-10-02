# Scoring support from partial RNA reads

`--partial-read-support` adds a separate, catalog-relative support score to
protein-hypothesis exports. It reconsiders original reads at the variant,
including partial reads that were not used in a particular assembly. The
existing reconstruction, Isovar ordering, raw RNA counts and cell compatibility
assignments remain available unchanged.

```sh
isovar protein-hypotheses --vcf variants.vcf --bam rna.bam \
    --sample-id sample-T1 --partial-read-support --output proteins.json
```

For pooled reconstruction and per-cell scores, also add
`--cell-reconstruction sample` (or `library`, `cell`, `group`). The same
comparison settings apply to both outputs:

```sh
isovar protein-hypotheses --vcf variants.vcf --bam rna.bam \
    --sample-id sample-T1 --cell-reconstruction sample \
    --partial-read-support --support-max-edits 1 --support-min-overlap 10 \
    --output proteins.json
```

The defaults are one anchored edit and ten compared bases. Existing
`--cell-max-edits` and `--cell-min-overlap` names remain aliases. These are
explicit compatibility settings, not calibrated platform-specific error rates.

## How the score works

Every scored fragment has one unit of support. If it is compatible with K
hypotheses, each receives 1/K. The score of a hypothesis is the sum of those
contributions. A fragment means a source-scoped read-group/QNAME pair: paired
mates, merged mates, repeated records and alternate placements cannot produce
extra units. These units are not deduplicated molecules or UMIs.

| Original evidence | Contribution to A | Contribution to B |
| --- | ---: | ---: |
| A short read distinguishing A under the comparison policy | 1 | 0 |
| A read observing only sequence shared by A and B | 0.5 | 0.5 |
| A reference-allele read | 0 | 0 |

Thus nine shared fragments and one A-specific partial fragment give A a score
of 5.5 and B a score of 4.5. Extending candidates outside the observed window
does not divide the read's contribution by the new candidate length. All
candidates within the configured tolerance remain alternatives: an exact match
and a one-edit match remain ambiguous when both pass that policy.

The matcher uses the observed prefix/allele/suffix, the focal reference-allele
competitor and compatible transcript paths. The recorded comparison policy is
`anchored_edit_compatibility.v2`: flanks align outward from the focal allele
until either sequence ends, allowing free distal overhangs while charging
internal and focal-adjacent edits. A fragment's jointly compatible
RNA candidates are the intersection across its observed mates. Alternative
placements must agree on their compatibility sets. Conflicting observations,
disagreeing placements, or any observation lacking a usable comparison leave
the fragment unscored, with its observations and reason retained.

RNA and protein scores are computed separately. For protein scores, the joint
RNA candidate set is projected to distinct protein hypothesis IDs **before**
dividing the unit. Two synonymous RNAs coding for the same protein therefore
do not double its support or dilute its share against another protein.

Scores depend on the supplied candidate catalog. `catalog_complete` says
whether all discovered hypotheses were retained under the recorded cap; it does
not establish that every possible expressed sequence was discovered. Local
windows and longer hypotheses remain visible, including the existing
`covered_by` relationships. Adding alternatives can change allocations.

## JSON and TSV

Each event gains `partial_read_support` with:

- The comparison/allocation policy and catalog completeness.
- Original `fragments`, with scoped hashed fragment IDs, RNA evidence-set IDs,
  comparisons, compatible RNA/protein IDs and assigned weights. Legacy
  QNAME-only observations retain null exported identities.
- `input_fragments`, `scored_fragments` and `unscored_fragments`.
- `candidate_scores` and `protein_scores`, each with `fractional_fragments`,
  `unique_fragments`, `ambiguous_fragments` and a dense `rank` (ties share rank).

Protein rows also carry their `partial_read_support` summary. The opt-in TSV
adds `protein_fractional_fragments`, `protein_unique_support_fragments`,
`protein_ambiguous_support_fragments` and `protein_support_score_rank`. The
existing `isovar_rank` and row order remain unchanged. Summed RNA scores and,
separately, summed protein scores equal the number of scored fragments, up to
floating-point export precision. Do not add the two resolutions together.

Cell reconstruction uses the same scorer within eligible biological sample
boundaries. Each resolved cell gains a score summary; individual fragment
assignments link to its scoped cell ID. Fragments without consistent cell
labels remain visible through `unresolved_cell_fragments`. A cell with one
partial read can contribute a full or fractional unit even when discovering
the candidate required reads pooled from many cells.

## Quality and interpretation

Adding `--min-base-quality N` annotates scored fragments as `passed`, `failed`
or `unassessed` under the existing [focal-allele quality policy](base-quality.md).
Any failed source read fails the fragment category; otherwise any unassessed
source read leaves it unassessed. Summaries include
`base_quality_fractional_fragments` using the same allocations. Missing QUAL
leaves the score unchanged and never becomes a fabricated Q0 or passing score.
These categories assess the focal allele (anchors for empty alleles), not the
flanking features that distinguish whole RNA/ORF candidates.

This is a transparent equal-allocation evidence score. It is not the
abundance/error model and posterior inference used by tools such as
[RSEM](https://pmc.ncbi.nlm.nih.gov/articles/PMC3163565/), nor an expression
probability, cell prevalence estimate or automatic pooling policy. Calibration
and quality-weighted likelihoods require additional modeling and controls.

Since 1.42.1 the shared anchored matcher corrects the equal-length clipping bug
in [#428](https://github.com/openvax/isovar/issues/428), which could count one
flanking insertion as two edits. Each comparison records `flank_alignments`
with all optimal `[read_bases, candidate_bases]` endpoint pairs, measured outward
from the allele, plus `compared_read_bases`. Coverage uses the minimum span
across tied endpoints, while transcript-path and A/C/G/T checks use the maximum
span. This prevents an endpoint tie from claiming certain full-candidate
coverage or hiding incompatible evidence. The focal allele remains globally
aligned, and unscored fragments stay available for inspection.

## Python API

```python
from isovar import export_protein_hypotheses, reconstruct_cell_groups

export = export_protein_hypotheses(
    results, sample_id="sample-T1", source="rna.bam",
    partial_read_support=True, support_max_edits=1, support_min_overlap=10,
)

cells = reconstruct_cell_groups(
    results, alignment_file, sample_id="sample-T1", source="rna.bam",
    mode="sample", partial_read_support=True, max_edits=1, min_overlap=10,
)
```

Validation includes explicit informative/shared-read examples, conserved
fragment budgets, mates and alternate placements, synonymous projection,
missing qualities, both transcript strands, sample boundaries and default
export parity. Independent public Illumina, ONT and PacBio allele ledgers and
original ONT cell reads with pinned transcript references run offline.
