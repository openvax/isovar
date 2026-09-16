# Mutation-evidence figures

Install the optional renderer (existing Isovar consumers need no new dependency):

```sh
pip install 'isovar[plot]'
```

Plot one mutation from your BAM and reference annotation:

```sh
isovar-plot --variant 5 141526090 C A --genome GRCh38 --bam tumor.bam \
  --compare-assembly --output-dir figures

isovar-plot --variant 5 141526090 C A --genome GRCh38 --bam tumor.bam \
  --view protein --compare-assembly --output-dir figures

isovar-plot --variant 5 141526090 C A --genome GRCh38 --bam tumor.bam \
  --view transcripts --output-dir figures

isovar-plot --variant 5 141526090 C A --genome GRCh38 --bam tumor.bam \
  --view assembly --output-dir figures
```

The default `--view all` combines protein, coverage, read-overlap and local
transcript panels. Existing variant-file, read-collection and translation
options work here too. Select exactly one mutation per invocation.
`--compare-assembly` changes **only** overlap assembly between the two runs;
mate merging and every other setting remain the same. This is not a comparison
of raw single-end reads against paired-end reads.

Each invocation creates a new directory; it never overwrites an earlier run:

```text
figures/
  YYYY-MM-DD_HH-MM-SS-microsecondsZ/     # UTC
    DIAPH1-5-141526090-C-A/
      overview.svg                    # vector graphics with editable text
      overview.png                    # opaque white, 300 dpi by default
      evidence.json                   # sequences, coordinates, support, settings
```

Use `--dpi` to change PNG resolution and `--max-rows` to limit displayed span
groups/transcript models. Neither reduces the reads used for analysis. A row
marked ×N represents N post-merge read objects with the same clipped span.
Coverage uses all supporting objects. Read names are not exported.

## Reading the figure

- Protein rows are aligned at the mutation, N terminus to C terminus. Orange
  marks altered residues or an in-frame deletion boundary. A frameshift marks
  the altered downstream region. Contexts over 80 aa use bars; complete
  sequences are always in `evidence.json`.
- RNA rows use transcript-oriented cDNA offsets: upstream is left even for a
  minus-strand transcript. Coverage and overlap panels share the same scale.
- One actual supporting cDNA witness is shown for each selected protein;
  protein support can include other translations of the same protein. The
  overlap panel shows the first mode with a protein, normally assembly on.
  Witness support and overall protein support are recorded separately.
- Transcript models use **forward genomic coordinates**, with strand arrows.
  Exons are clipped to the local observed region; long intronic gaps are
  compressed and marked `//`. Junction counts come from retained CIGAR N
  evidence for the displayed witness, not from assuming the annotation is true.
- The figure title uses Varcode's normalized 1-based allele representation.
  Genomic track boundaries are 0-based, half-open. An anchored VCF deletion
  can therefore have a different displayed start after its anchor is removed.

Models that contribute the same protein do not uniquely identify an isoform.
The figure is a rendering of Isovar's output, not an independent validation of
its sequence, full-length transcript, expression, peptide presentation or
clinical suitability. No extra sequence is filled in from the reference.
Like `isovar-protein-sequences`, it shows candidates before `run_isovar`'s
result-level filters; producing a plot does not mean a variant passed them.

Known interpretation limitations remain tracked in
[#264](https://github.com/openvax/isovar/issues/264) (secondary alignments in
mate merging) and [#265](https://github.com/openvax/isovar/issues/265) (upstream
indels when assigning the reading frame). Plotting does not repair them.

## Reproduce the osteosarc examples

Browse the [recorded figure gallery](../figures/osteosarc/2026-09-16_01-26-26-429446Z/README.md).

From a repository checkout with the plotting extra installed:

```sh
python -m examples.osteosarc_assembly_figures --output-dir figures/osteosarc
```

This uses checksum-pinned, original primary-only RNA fixtures and offline
Ensembl 87 models. No download, altered read, synthetic read or realignment
is needed. It writes a UTC-stamped batch containing three variant directories,
individual captions and an overall manifest. Both modes must pass independent
translation and transcript-reference checks before figures are written.

| Example | Assembly off → on | Mutation-overlapping 25-mers |
| --- | --- | --- |
| DIAPH1 | 40 → 49 aa | 16 → 25 |
| SLC25A12 | 44 → 49 aa | 20 → 25 |
| TECPR1 | 37 → 49 aa | 13 → 25 |

In each example, no single read object spans the full assembled cDNA witness.
Assembly connects overlapping observations to recover missing context. The
shorter no-assembly sequence is correct, but incomplete for those additional
peptide windows. Selected fixture counts are not full-sample abundance estimates.

The public data are CC0; see the
[source and selection documentation](../tests/data/osteosarc/README.md).
The batch records the exact source BAM URLs, original alleles, fixture/reference
checksums, effective settings and independent checks, without embedding read names.

## Python

```python
from isovar.visualization import (
    collect_visualization_data, save_variant_figures, timestamped_run_directory,
)

# variant and read_evidence are ordinary Isovar pipeline inputs.
data = collect_visualization_data(variant, read_evidence, compare_assembly=True)
run = timestamped_run_directory("figures")
directory = save_variant_figures(data, run)
```

`collect_visualization_data` does not import Matplotlib. Its `creator_kwargs`
argument accepts existing `ProteinSequenceCreator` settings. For further
styling, `plot_variant_evidence(data)` returns a Matplotlib Figure without
selecting a GUI backend or changing global plotting defaults.
