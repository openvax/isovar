# Reading and filtering results

This page explains what `run_isovar` returns, what each result and table column
means, and how filters mark results. Read it after the
[quick start](../README.md#quick-start), before relying on Isovar's output.

## Reading the results

`run_isovar` returns one `IsovarResult` per input variant, in input order. Each
holds the RNA evidence at that variant's locus and any mutant protein sequences
assembled for it.

| `IsovarResult` property | Table column | Meaning |
|---|---|---|
| `variant` | `variant`, `chr`, `pos`, `ref`, `alt` | The input variant |
| `top_protein_sequence` | `protein_sequence` and other `protein_sequence_*` columns | The top-ranked mutant protein assembled from RNA, or none |
| `sorted_protein_sequences` | | Every protein kept, ranked; one by default |
| `num_alt_reads`, `num_alt_fragments`, `fraction_ref_reads`, ... | Same names | Read evidence for the ref, alt and other alleles |
| `predicted_effect` | `predicted_effect` | Varcode's selected comparison with its retained alternatives and supplied patient-baseline evidence |
| `has_protein_changing_predicted_effect` | Same name | Whether any retained candidate predicts a protein change; unknown if unresolved and none is known to change |
| `passes_all_filters`, `filter_values` | `passes_all_filters`, `filter:<name>` | Filter outcomes |

By default Isovar aims for 49 residues of protein: 2K−1 for peptide size
K = 25, enough to contain every 25-mer overlapping a single changed residue.
The actual length depends on RNA coverage.
[Protein context selection](protein-selection.md)
explains the target and how the top protein is chosen.

A read is one sequenced segment, such as one mate of a pair; a fragment is one
template, so both mates of a pair count once. Fragment counts are usually the
ones to filter on. Neither is a molecule count; see
[collecting and counting RNA reads](read-collection.md).

## Filtering results

`run_isovar` evaluates filters on each result; a failing result is kept, with
`False` in its `filter_values` dictionary and in `passes_all_filters`. In the
DataFrame and CSV each filter becomes a `filter:<name>` column. Filters never
remove a result or a row.

`filter_thresholds` maps names like `'min_num_alt_reads'` or
`'max_fraction_other_fragments'` to numbers. The text after `min_` or `max_`
names a numeric property of `IsovarResult`, and most read-evidence properties
follow the pattern `{num|fraction}_{ref|alt|other}_{reads|fragments}`. For
example, this requires at least 10 alt reads and at most 25% of fragments
supporting other alleles:

```python
from isovar import run_isovar

isovar_results = run_isovar(
    variants="cancer-mutations.vcf",
    alignment_file="tumor-rna.bam",
    filter_thresholds={"min_num_alt_reads": 10, "max_fraction_other_fragments": 0.25})

for isovar_result in isovar_results:
    print(isovar_result.variant, isovar_result.passes_all_filters)
```

`filter_flags` names boolean properties of `IsovarResult`; prefix one with
`not_` to negate it, as in `not_protein_sequence_matches_predicted_mutation_effect`.

Omitting either argument applies the defaults in
[`default_parameters.py`](../isovar/default_parameters.py):
`DEFAULT_FILTER_THRESHOLDS`, and `DEFAULT_FILTER_FLAGS`, which are
`has_protein_changing_predicted_effect`,
`has_mutant_protein_sequence_from_rna` and `protein_sequence_contains_mutation`.
Passing a value replaces the corresponding defaults; to change one threshold,
copy `DEFAULT_FILTER_THRESHOLDS` and update it.

`predicted_effect_modifies_protein_sequence` continues to describe the selected
comparison alone. The default candidate-aware filter also considers alternatives
when that preferred result is intronic or synonymous. Unknown change status
fails this positive filter. See [prediction comparisons](reference-annotation.md#prediction-comparisons).
