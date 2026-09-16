# Mutation-evidence figures

Install the optional renderer (existing Isovar consumers need no new dependency):

```sh
pip install 'isovar[plot]'
```

Plot one mutation from your BAM and reference annotation:

```sh
isovar plot --variant 5 141526090 C A --genome GRCh38 --bam tumor.bam \
  --compare-assembly --output-dir figures

isovar plot --variant 5 141526090 C A --genome GRCh38 --bam tumor.bam \
  --view protein --compare-assembly --output-dir figures

isovar plot --variant 5 141526090 C A --genome GRCh38 --bam tumor.bam \
  --view transcripts --output-dir figures

isovar plot --variant 5 141526090 C A --genome GRCh38 --bam tumor.bam \
  --view reads --output-dir figures
```

The default `--view all` exports **four individual figures** (protein, coverage,
read overlaps and transcripts) plus a combined overview. Individual panels use
a wider 16-inch canvas with editable SVG and a 600-dpi PNG (9,600 pixels wide).
The overview PNG is a 300-dpi preview (4,800 pixels wide); its SVG and the PDF
remain fully vector. Brief
support notes sit in a separate right margin; definitions and detailed
provenance live in `caption.md` and `evidence.json`, not over the data.

Select `--view protein`, `coverage`, `reads` or `transcripts` to export only one.
The existing `--view assembly` still produces combined coverage/read panels,
now with their separate figures too. The `isovar-plot` alias remains supported.
Existing variant-file, read-collection and translation
options work here too. Select exactly one mutation per invocation.
`--compare-assembly` changes **only** overlap assembly between the two runs;
mate merging and every other setting remain the same. This is not a comparison
of raw single-end reads against paired-end reads.

Each invocation creates a new directory; it never overwrites an earlier run:

```text
figures/
  YYYY-MM-DD_HH-MM-SS-microsecondsZ/     # UTC
    DIAPH1-5-141526090-C-A/
      overview.{svg,png}               # combined overview
      protein.{svg,png}                # standalone, large protein panel
      coverage.{svg,png}               # standalone coverage panel
      reads.{svg,png}                  # standalone read-overlap panel
      transcripts.{svg,png}            # standalone transcript panel
      caption.md                      # definitions and interpretation
      all-figures.pdf                  # one standalone panel per page
      evidence.json                   # sequences, coordinates, support, settings
```

All backgrounds are opaque white. Use `--dpi` to change PNG resolution (overview
previews are capped at 300 dpi) and `--max-rows` to limit displayed span
groups/transcript models. Neither reduces the reads used for analysis. A row
marked ×N represents N post-merge read objects with the same clipped span.
Coverage uses all supporting objects. Read names are not exported.

## Reading the figure

- Protein rows are aligned at the mutation, N terminus to C terminus. Orange
  marks altered residues or an in-frame deletion boundary. A frameshift marks
  the altered downstream region. Contexts over 80 aa use bars; complete
  sequences are always in `evidence.json`.
- Varcode rows predict **reference transcript + nominated edit only**, on the
  same annotation. They are not RNA reconstructions. Equal local predictions
  share a row, with every transcript and effect label mapped in `caption.md`.
  The alignment origin is the genomic edit's codon, not Varcode's potentially
  repeat-normalized effect-label boundary. Missing Isovar context stays missing.
  Black transcript boxes contribute to the displayed RNA protein; gray boxes
  are prediction-only models. Alignment outside the displayed cDNA can affect
  compatibility, so this local track does not explain every model exclusion.
- RNA rows use transcript-oriented cDNA offsets: upstream is left even for a
  minus-strand transcript. Coverage and overlap panels share the same scale.
- **Reconstructed cDNA** means one actual RNA-derived sequence that produces
  the selected protein. Several different reconstructions can yield the same
  protein; this track does not combine their reads or establish a unique isoform.
  Protein support can include other translations of the same protein. The
  overlap panel shows the first mode with a protein, normally assembly on.
  Reconstruction support and overall protein support are recorded separately.
  The JSON key `witness` is retained for compatibility with existing consumers.
- Transcript models use **forward genomic coordinates**, with strand arrows.
  ENST labels include transcript names from the same annotation when available.
  Angled gray connectors join adjacent exon boundaries; black connectors on the
  separately labeled RNA junction row show observed splice evidence.
  Exons are clipped to the local observed region; long genomic gaps outside
  the displayed exon union are compressed and marked `//`. These can include
  unannotated flanks, not just introns. Junction counts come from retained CIGAR N
  evidence for the displayed reconstruction, not from assuming the annotation is true.
- The figure title uses Varcode's normalized 1-based allele representation.
  Genomic track boundaries are 0-based, half-open. An anchored VCF deletion
  can therefore have a different displayed start after its anchor is removed.

Models that contribute the same protein do not uniquely identify an isoform.
The figure is a rendering of Isovar's output, not an independent validation of
its sequence, full-length transcript, expression, peptide presentation or
clinical suitability. No extra sequence is filled in from the reference.
Coordinate conversions follow the [Ensembl GTF definition](https://www.ensembl.org/info/website/upload/gff.html);
observed introns use CIGAR N per the [SAM specification](https://samtools.github.io/hts-specs/SAMv1.pdf).
Like `isovar protein-sequences`, it shows candidates before `run_isovar`'s
result-level filters; producing a plot does not mean a variant passed them.

Known interpretation limitations remain tracked in
[#264](https://github.com/openvax/isovar/issues/264) (secondary alignments in
mate merging) and [#265](https://github.com/openvax/isovar/issues/265) (upstream
indels when assigning the reading frame). Plotting does not repair them.

## Reproduce the osteosarc examples

Browse the [expanded figure gallery](../figures/osteosarc/2026-09-16_05-43-34-415330Z/README.md)
or the [original three-example, 12-page vector PDF](../figures/osteosarc/2026-09-16_05-23-24-826750Z/three-examples-12-pages.pdf).

From a repository checkout with the plotting extra installed:

```sh
python -m examples.osteosarc_assembly_figures --output-dir figures/osteosarc
```

This uses checksum-pinned, original primary-only RNA fixtures and offline
Ensembl 87 models. No download, altered read, synthetic read or realignment
is needed. It writes eleven example directories, individual captions and a
manifest. Every returned RNA protein must pass independent translation/frame/
interval checks; every displayed Varcode prediction must match the independent
reference-plus-edit translation. RNA equality to that prediction is a separate
check, not an assumption.

| Example | Assembly off → on | Mutation-overlapping 25-mers |
| --- | --- | --- |
| DIAPH1 | 40 → 49 aa | 16 → 25 |
| SLC25A12 | 44 → 49 aa | 20 → 25 |
| TECPR1 | 37 → 49 aa | 13 → 25 |

In these first three examples, no single read object spans the full reconstructed cDNA.
Assembly connects overlapping observations to recover missing context. The
shorter no-assembly sequence is correct, but incomplete for those additional
peptide windows. Selected fixture counts are not full-sample abundance estimates.

Additional examples:

- **DYNC1H1:** both the V314I and Q3267H loci, kept separate.
- **H1-2:** 15-nt in-frame deletion; the pinned annotation calls the gene
  HIST1H1C. Varcode repeat-normalizes its protein deletion label differently
  from the source, but its predicted sequence is checked independently.
- **MAP2:** 22-nt frameshift deletion, one alternate RNA object, no protein
  under default coverage requirements. Two reference predictions are shown
  without fabricating RNA support.
- **NAV2:** two transcripts contribute the current 47-aa RNA sequence. It ends
  at W, before the reference-only tracks diverge into LR versus VN. Eight
  annotation models are retained, seven with protein predictions. The retained
  cDNA has no junction; the local display does not establish one isoform or
  justify adding either predicted terminal extension. This supersedes the
  earlier pre-splice-restriction 49-aa V/L comparison.
- **NTF3:** Varcode's nominated A>G predicts K56R (K69R in the other model),
  but RNA encodes S. The source-linked adjacent change produces a compound
  AG>GT allele. Independent translation checks pass, while the single-edit
  reference comparison fails as expected. This is not a claim that the second
  change is germline, nor that assembly is required (both modes recover S).
- **PIP5K1A T1 ONT versus T1 Illumina:** complete regional inputs, not selected
  fixtures. ONT gives 46 aa in both modes; Illumina has one alternate object
  and no protein at the two-read floor. This is a library-level support
  difference, **not a controlled read-length comparison**. See the
  [pinned original-read comparison and reproduction recipe](../tests/data/osteosarc/figure_comparisons/README.md).

### DIAPH1 junction-track walkthrough

The five named transcript models share the displayed intron. Gray angled lines
join annotated exon boundaries. The separate black RNA junction connects
`[141524229, 141526037)` on the forward, 0-based genomic axis; its label **10**
counts retained post-merge observations, not intron length or independent
molecules. DIAPH1 is on the minus strand, so its arrows point left, while the
cDNA/protein panels still run 5′→3′ / N→C left-to-right. `//` means the genomic
intron is compressed for display. The shared junction supports a local splice
path but does not choose one of the five isoforms.

### Structural variants

This gallery does not establish reconstruction of a real osteosarc SV. The
pinned 44-allele cohort contains SNVs and small deletions, not a validated
breakpoint/fusion example. Literal-sequence indels and symbolic SVs are different
inputs: Isovar has no general breakpoint/END/orientation-aware RNA reconstruction
path. In Varcode 7.0.0 and 9.0.0, a symbolic `StructuralVariant` exposes placeholder
bases that Isovar's small-allele helpers currently interpret incorrectly;
see [#270](https://github.com/openvax/isovar/issues/270). Do not interpret a
symbolic SV result as a validated rearranged protein. Fixing that integration
and validating real breakpoint-spanning RNA are separate work, not hidden
inside the renderer.

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
control, `include_reference_predictions=False` omits Varcode tracks without
changing RNA analysis. The JSON additions are `reference_predictions`,
`varcode_version` and per-transcript `rna_supported` (contributes to a displayed
protein, not unique isoform identification). For further
styling, `plot_variant_evidence(data)` returns a Matplotlib Figure without
selecting a GUI backend or changing global plotting defaults.
