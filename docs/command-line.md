# Command line

This page lists every `isovar` command, the conventions its output tables
share, and how to change options, logging and error handling. The pipeline
command is:

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

### Pipeline

| Command | Output |
|---|---|
| `isovar run` | One row per variant: read evidence, top protein sequence, predicted effect and filters |

### Table exports

| Command | Output |
|---|---|
| `isovar protein-sequences` | Ranked candidate protein sequences (`--max-protein-sequences-per-variant 0` keeps all) |
| `isovar translations` | Every translation of each assembled cDNA in each compatible reading frame, before grouping |
| `isovar variant-sequences` | Assembled cDNA sequences supporting each variant |
| `isovar variant-reads` | Reads supporting each variant's alt allele |
| `isovar allele-reads` | All reads overlapping each variant |
| `isovar allele-counts` | Read and fragment counts for the ref, alt and other alleles; with `--cell-umi-labels`, cells and cell/UMI labels too |
| `isovar reference-contexts` | Reference sequence and reading frame around each variant (no BAM needed) |

### RNA hypotheses and interpretation

| Command | Output |
|---|---|
| `isovar protein-hypotheses` | Every protein hypothesis, its translations and scoped RNA evidence, as JSON and TSV ([guide](protein-hypotheses.md)) |
| `isovar variant-rna-paths` | Original ALT/REF RNA paths assembled before coding-model assessment, as JSON ([guide](variant-rna-paths.md)) |
| `isovar allele-interpretations` | Which competing representations of one locus the RNA reads carry, as JSON ([guide](allele-interpretations.md)) |
| `isovar rna-features` | Read/fragment counts and optional library/feature normalization for explicit RNA paths ([guide](https://github.com/openvax/isovar/blob/master/RNA_FEATURE_SUPPORT.md)) |
| `isovar rna-origin` | Nuclear/mitochondrial context evidence against an explicit panel, as JSON ([guide](rna-origin.md)) |

### Structural variants and fusions

| Command | Output |
|---|---|
| `isovar sv-source` | Validate original SV definitions, caller claims and record-specific PCR evidence ([guide](sv-source-validation.md)) |
| `isovar sv-rna` | Exploratory RNA paths around one nominated SV, as JSON ([guide](sv-rna.md)) |
| `isovar sv-audit` | Reconstruct and reconcile a source-separated SV inventory ([guide](sv-rna-audit.md)) |
| `isovar fusion-rna` | Validated junction evidence and frames for a supplied fusion, as JSON (compatibility alias: `fusion`; [guide](fusion.md)) |

### Viral and isoform evidence

| Command | Output |
|---|---|
| `isovar viral-rna` | Observed viral RNA recovered from combined-reference alignments without a GTF (compatibility alias: `viral`; [guide](viral-sequences.md)) |
| `isovar isoforms` | Targeted antigen (for example CTA) isoforms reconstructed and quantified (`reconstruct`, `quantify`; [guide](antigen-isoforms.md)) |

### Proteome and patient background

| Command | Output |
|---|---|
| `isovar rna-proteome` | Recorded RNA hypotheses adapted to an ORFdata Proteome and FASTA (compatibility alias: `proteome`; [guide](rna-proteome.md)) |
| `isovar patient-background-pairs` | Mutant RNA proteins paired with declared patient-background baselines ([guide](patient-background-pairs.md)) |
| `isovar normal-evidence` | Normal genotype/reference-confidence states and background alleles from a reviewed source inventory ([guide](matched-normal-evidence.md)) |

### Reports and figures

| Command | Output |
|---|---|
| `isovar plot` | Protein, coverage, read-overlap and transcript figures for one mutation ([guide](visualization.md)) |
| `isovar evidence-report` | Visual PDF and per-call TSV for a variant VCF plus indexed tumor RNA |
| `isovar antigen-report` | Source-separated PDF and JSON from small-variant, SV, fusion, viral and antigen-isoform exports ([guide](antigen-evidence-report.md)) |
| `isovar cohort-report` | Reproducible, source-separated reconstruction and per-sample PDFs from a pinned manifest (compatibility alias: `cohort`; [guide](cohort-reports.md)) |

### Calibration

| Command | Output |
|---|---|
| `isovar calibrate-phasing` | Phasing error bounds calibrated from truth-labeled control BAMs ([guide](phasing.md)) |

All subcommands are invoked through the single `isovar` executable (`isovar COMMAND ...`).

## Output tables

`isovar evidence-report` writes a PDF and a tab-separated per-call file. Reports
are included in the normal installation:

```sh
pip install isovar
isovar evidence-report --vcf somatic-variants.vcf.gz --bam tumor-rna.bam \
    --sample-id tumor-1 --output evidence.pdf
```

Use `--nc-orf-catalog` to add reference ncORF coding-block candidates and
`--self-fasta` to rank mutation-overlapping 9-mers absent from a reference
protein FASTA. The ncORF catalogue is checked against the variant reference;
the supplied FASTA path and SHA-256 are recorded so its reference identity can
be reviewed. The
command reports small variants; nominated SV breakends need a separate SV RNA
screen. The PDF distinguishes source DNA predictions, original RNA allele
reads and local Isovar protein windows. It does not infer full-length isoforms.
It includes exon-involved ALT-read distributions with and without zero-support
calls. Use `isovar antigen-report` to compose other source-specific exports,
including separate viral fragment distributions.

The following table conventions apply to the other table commands (`run` and the
table exports).

Every table starts with the same variant key: `variant` (as in
`chr9 g.82927102G>T`) and `chr`, `pos`, `ref` and `alt` as given in the input,
so tables from different commands can be joined.

- Default output filenames follow `isovar-<command>.csv` (and `isovar-results.csv` for `isovar run`).
- Output paths ending in `.tsv` or `.tab` (optionally followed by `.gz`/`.bz2`/`.xz`) write tab-separated tables; `--output -` writes to stdout while logs remain on stderr.
- Gene lists use `gene_names`, transcript lists use `transcript_ids` / `transcript_names`, and candidate support uses `num_supporting_reads` / `num_supporting_fragments`.
- `allele-reads` and `variant-reads` include `allele_class` (`ref`, `alt`, or `other`), `mate` (`R1`, `R2`, `R1+R2`, or `unpaired`), and `is_secondary` so rows reconcile with `allele-counts`.
- Counts are `num_*` columns.
- Lists are `;`-separated.
- Numbers are written with six significant digits.
- An empty result still has its header.

The CLI applies the same default filters as `run_isovar`. Filters only set the
`filter:*` columns and `passes_all_filters`; they never remove rows. See
[filtering results](results.md#filtering-results).

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

## CLI conventions

Isovar enforces a consistent command and flag vocabulary across every subcommand:

- **Command naming:** Data-producing commands use lowercase hyphenated artifact/domain noun phrases (plural nouns for table exports such as `protein-sequences` and `allele-counts`; domain-artifact pairs such as `sv-rna`, `fusion-rna`, `viral-rna`, `rna-proteome`, `evidence-report`, `antigen-report`, and `cohort-report`). Imperative verbs are reserved for top-level actions (`run`, `plot`, `calibrate-phasing`). Deprecated bare-noun aliases (`fusion`, `viral`, `proteome`, `cohort`) remain accepted with a one-line warning on stderr.
- **Flat vs. nested subcommands:** Commands are flat (`isovar COMMAND ...`) by default. Nested subcommands are used only for multi-step workflows that operate sequentially on the same domain artifact (`isovar isoforms reconstruct` followed by `isovar isoforms quantify`).
- **Output flags:** `-o`/`--output FILE` always denotes a single output file, and `--output-dir DIR` always denotes a directory. Companion files use `--output-tsv`, `--output-json`, or `--output-fasta`. Legacy spellings (`--plot-dir`, `--summary-json`, `--evidence-json`, `--per-call`, `--fasta`) remain accepted as hidden compatibility aliases.
- **Input vocabulary:** Shared inputs use one spelling across commands: `--bam` (indexed RNA BAM/CRAM), `--vcf` (variants), `--input` (structured JSON/upstream export), `--manifest` (multi-source/cohort manifest JSON), `--genome` (Varcode/PyEnsembl reference), `--genome-fasta` (local FASTA; hidden alias `--reference-fasta`), `--gtf` (local transcript annotation), `--sample-id` (sample scope), `--source` (read-set/library identity; hidden alias `--source-id`), `--library-cache` (reusable RNA index directory; hidden alias `--library-cache-dir`), and `--scratch-dir` (temporary disk-index directory; hidden alias `--sv-scratch-dir`).
- **Short flags:** Only globally unambiguous single-letter flags are defined: `-h`/`--help`, `-o`/`--output`, `-f`/`--force` (overwrite existing outputs on overwrite-guarded commands), and `-q`/`--quiet` (suppress INFO progress logs on stderr).
- **Alignment with sibling OpenVax tools:** Isovar adopts `-o`/`--output`, `-f`/`--force`, `-q`/`--quiet` alongside `--log-level`, explicit `--version` help text, and JSON output across all structured-analysis commands. Deliberate differences: Isovar keeps `--genome` (via Varcode, which binds species to an explicit Ensembl release or assembly) instead of a separate `--species` flag, and leaves reference-catalog management (`list`, `status`, `download`, `install`, `--json`) in `orfdata` and `pyensembl` rather than duplicating those commands in `isovar`.

## Logging and errors

Progress messages go to stderr. `--log-level DEBUG` adds per-candidate detail,
and `-q`/`--quiet` (or `--log-level WARNING`) keeps runs quiet.

Unusable inputs and out-of-range options stop with a one-line error and exit
status 2. Examples are a missing file, an unindexed BAM, SAM instead of BAM, a
missing output directory and malformed JSON.

## Shell tab completion

`isovar` supports tab completion for subcommands and options via `argcomplete`.
Turn it on once for every OpenVax command-line tool:

```sh
activate-global-python-argcomplete
```

Or enable it per tool:

```sh
eval "$(register-python-argcomplete isovar)"
```
