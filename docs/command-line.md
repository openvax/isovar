# Command line

```sh
isovar run \
    --vcf somatic-variants.vcf \
    --bam rnaseq.bam \
    --output isovar-results.csv
```

`isovar --help` lists the subcommands; each subcommand's `--help` lists its
options and defaults, which match the Python API. `isovar --vcf ... --bam ...`
(without a subcommand) also runs the pipeline, and `python -m isovar` works too.
`--genome` chooses the annotation; see
[choosing the reference annotation](reference-annotation.md).

## Commands

| Command | Output |
|---|---|
| `isovar run` | One row per variant: read evidence, top protein sequence, predicted effect and filters |
| `isovar protein-sequences` | Ranked candidate protein sequences (`--max-protein-sequences-per-variant 0` keeps all) |
| `isovar protein-hypotheses` | Every protein hypothesis, its translations and scoped RNA evidence, as JSON and TSV ([guide](protein-hypotheses.md)) |
| `isovar translations` | Every translation of each assembled cDNA in each compatible reading frame, before grouping |
| `isovar variant-sequences` | Assembled cDNA sequences supporting each variant |
| `isovar reference-contexts` | Reference sequence and reading frame around each variant (no BAM needed) |
| `isovar allele-counts` | Read and fragment counts for the ref, alt and other alleles; with `--cell-umi-labels`, cells and cell/UMI labels too |
| `isovar allele-reads` | All reads overlapping each variant |
| `isovar variant-reads` | Reads supporting each variant's alt allele |
| `isovar plot` | Protein, coverage, read-overlap and transcript figures for one mutation ([guide](visualization.md)) |
| `isovar allele-interpretations` | Which competing representations of one locus the RNA reads carry, as JSON ([guide](allele-interpretations.md)) |
| `isovar sv-rna` | Exploratory RNA paths around one nominated SV, as JSON ([guide](sv-rna.md)) |
| `isovar fusion` | Validated junction evidence and frames for a supplied fusion, as JSON ([guide](fusion.md)) |
| `isovar audit-bundle` | Portable audit inputs, verified backups and offline reconstruction replay ([guide](audit-bundles.md)) |
| `isovar uorf-evidence` | Query and export the offline human uORF evidence reference ([guide](uorf-evidence.md)) |

Except `isovar run`, the table and plot commands also install as hyphenated
scripts such as `isovar-protein-sequences` and `isovar-plot`.

## Output tables

Every CSV starts with the same variant key: `variant` (as in
`chr9 g.82927102G>T`) and `chr`, `pos`, `ref` and `alt` as given in the input,
so tables from different commands can be joined.

- Counts are `num_*` columns.
- Lists are `;`-separated.
- Numbers are written with six significant digits.
- An empty result still has its header.

The CLI applies the same default filters as `run_isovar`. Filters only set the
`filter:*` columns and `passes_all_filters`; they never remove rows. See
[filtering results](../README.md#filtering-results).

## Changing options

For example, to use only primary alignments, include soft-clipped bases,
require at least three reads at every retained cDNA base, and decompress the BAM
with four threads:

```sh
isovar run --vcf somatic-variants.vcf --bam rnaseq.bam \
    --drop-secondary-alignments --use-soft-clipped-bases \
    --min-variant-sequence-coverage 3 --num-rna-decompression-threads 4 \
    --output isovar-results.csv
```

`--reference-context-size` applies only to `isovar reference-contexts`. The
protein commands size their reference context from the requested cDNA length
and minimum transcript prefix.

## Logging and errors

Progress messages go to stderr. `--log-level DEBUG` adds per-candidate detail,
and `--log-level WARNING` keeps runs quiet.

Unusable inputs and out-of-range options stop with a one-line error and exit
status 2. Examples are a missing file, an unindexed BAM, SAM instead of BAM, a
missing output directory and malformed JSON.
