# Mutation-evidence figures

`isovar plot` draws the evidence behind Isovar's result for one mutation:

- the RNA-derived protein, next to Varcode's reference prediction;
- RNA coverage along the reconstructed cDNA;
- the reads that overlap it; and
- the local transcript models, with the splice junctions the reads show.

## Quick start

```sh
pip install 'isovar[plot]'

isovar plot --variant 5 141526090 C A --genome GRCh38 --bam tumor.bam \
  --compare-assembly --output-dir figures
```

Each run writes a new UTC-stamped directory and never overwrites an earlier one:

```text
figures/
  YYYY-MM-DD_HH-MM-SS-microsecondsZ/
    DIAPH1-5-141526090-C-A/
      overview.{svg,png}      # all panels together
      protein.{svg,png}       # one file per panel
      coverage.{svg,png}
      reads.{svg,png}
      transcripts.{svg,png}
      all-figures.pdf         # one panel per page
      caption.md              # definitions and how to interpret the figures
      evidence.json           # sequences, coordinates, support and settings
```

Each panel is a 16-inch-wide figure, written as editable SVG and as a 600-dpi
PNG (9,600 pixels wide). The overview PNG is a 300-dpi preview; its SVG and the
PDF are fully vector. Backgrounds are opaque white. Brief support notes sit in
a right margin, while definitions and detailed provenance are in `caption.md`
and `evidence.json` rather than over the data. Read names are not exported.

## Options

| Option | Effect |
|---|---|
| `--view all` (default) | All four panels, plus the combined overview |
| `--view protein`, `coverage`, `reads` or `transcripts` | One panel |
| `--view assembly` | The combined coverage and read panels, plus their separate figures |
| `--compare-assembly` | Reconstruct twice, with overlap assembly off and on; nothing else changes |
| `--all-proteins` | Also write every recovered protein and frame alternative, paginated under `protein-alternatives/` |
| `--sample-label TEXT` | Name the sample and technology in the figures (default: the BAM filename) |
| `--dpi N` | PNG resolution (default 600; the overview preview is capped at 300) |
| `--max-rows N` | Show at most N read-span or transcript rows (default 20); analysis still uses every read |

Plot one mutation per run. The variant, read-collection and translation options
of the other commands work here too, and the command is also installed as
`isovar-plot`.

Some notes on these options:

- `--compare-assembly` changes **only** overlap assembly. The same reads are
  collected, and mate merging and every other setting stay the same. It does
  not compare single reads with paired reads. Sample and technology describe
  the input; assembly off and on are two reconstruction modes of that input.
- `--all-proteins` removes only the cap on protein results. Shorter contexts
  with identical sequence, frame and transcripts are grouped, without joining
  sequences or summing support. Ambiguous shorter contexts keep every compatible
  group. These are the results supported under the chosen thresholds, not every
  possible biological ORF.

## Reading the figure

**Protein panel.** Rows are aligned at the mutation, N terminus to C terminus.
Orange marks altered residues or an in-frame deletion boundary; for a
frameshift, the altered downstream region. Contexts longer than 80 aa are drawn
as bars, and complete sequences are always in `evidence.json`.

**Varcode rows** show the **reference transcript plus the nominated edit only**,
on the same annotation. They are predictions, not RNA reconstructions.
- Identical local predictions share a row; `caption.md` maps every transcript
  and effect label.
- Rows align at the codon of the genomic edit, not at Varcode's effect-label
  boundary, which can be shifted by repeat normalization.
- Black transcript boxes contribute to the displayed RNA protein; gray boxes are
  prediction-only models. Alignment outside the displayed cDNA can also affect
  compatibility, so this local view does not explain every excluded model.
- Missing Isovar context stays missing.

**RNA rows** use transcript-oriented cDNA offsets: upstream is on the left, even
for a minus-strand transcript. The coverage and read panels share one scale.
- **Reconstructed cDNA** is one actual RNA-derived sequence that produces the
  selected protein (the `witness` in `evidence.json`). Several reconstructions
  can give the same protein; this track neither combines their reads nor
  identifies a unique isoform.
- Protein support can include other translations of the same protein, so
  reconstruction support and protein support are recorded separately.
- The read panel shows the first mode that produced a protein, normally
  assembly on.
- A row marked ×N stands for N observations with the same clipped span, after
  mate merging and deduplication. Alternative alignments of one read cannot
  inflate support, and merged mates count once toward coverage. Distinct
  templates are not necessarily independent molecules.

**Transcript panel.** Models use **forward genomic coordinates**, with strand
arrows. Labels give the Ensembl transcript ID and, when the annotation has one,
the transcript name.
- Straight gray connectors join adjacent exon boundaries. Black connectors, on
  the separately labeled RNA junction row, show observed splicing. Junction
  counts come from the CIGAR `N` operations of the displayed reconstruction,
  not from assuming the annotation is true.
- Exons are clipped to the locally observed region. Long genomic gaps outside
  the displayed exons are compressed and marked `//`; these can be unannotated
  flanks as well as introns.

**Coordinates.** The title uses Varcode's normalized, 1-based allele. Genomic
track boundaries are 0-based and half-open, so an anchored VCF deletion can
show a different start once its anchor base is removed. Conversions follow the
[Ensembl GTF definition](https://www.ensembl.org/info/website/upload/gff.html),
and introns use CIGAR `N` as in the
[SAM specification](https://samtools.github.io/hts-specs/SAMv1.pdf).

### Example: the DIAPH1 junction track

The five named transcript models share the displayed intron, and gray straight
lines join their annotated exon boundaries. The black RNA junction connects
`[141524229, 141526037)` on the forward, 0-based genomic axis. Its label **10**
counts retained observations after mate merging, not intron length or
independent molecules. DIAPH1 is on the minus strand, so its arrows point left,
while the cDNA and protein panels still run 5′→3′ and N→C from left to right.
`//` marks the compressed intron. The shared junction supports a local splice
path but does not choose one of the five isoforms.

## How reads are handled

- **Missing base qualities.** A read without QUAL is kept. Alignment and
  sequence filters still apply, and its base confidence stays unknown. Where
  overlapping mates disagree and either quality is missing, the disagreement is
  not resolved by quality. `--require-base-qualities` discards such reads
  (`ReadCollector(use_reads_without_base_qualities=False)` in Python). Measured
  scores are never replaced; MAPQ and consensus read accuracy are not per-base
  Phred substitutes.
- **Soft clips.** `--use-soft-clipped-bases` keeps unaligned read ends. It does
  not realign an SV partner or recover an allele seen only in clips; CIGAR
  insertions are used either way. In a
  [matched osteosarc audit](../tests/data/osteosarc/figure_comparisons/SOFT_CLIPS.md),
  turning it on shortened some top protein contexts, because clipped sequence
  conflicted with aligned reads, so it stays off by default.
- **Reads ending at an insertion.** A read that ends at an insertion supports
  the reference allele only if both flanking reference bases are aligned. The
  absence of aligned inserted bases is not enough. Terminal insertions in the
  CIGAR are still used.
- **Filters.** Like `isovar protein-sequences`, the figures show candidates
  before `run_isovar`'s result filters. A plot does not mean the variant passed
  them.

## What a figure does not show

The figure renders Isovar's output. It is not an independent validation of the
sequence, a full-length transcript, expression, peptide presentation or
clinical suitability. Models that give the same protein do not identify one
isoform, and no sequence is filled in from the reference. Reading-frame
transfer and alternative-alignment handling are the production pipeline's,
including the fixes for [#264](https://github.com/openvax/isovar/issues/264)
and [#265](https://github.com/openvax/isovar/issues/265). A translation still
does not establish biological independence, a unique isoform or the absence of
sequencing errors.

## Osteosarc examples

Browse the [expanded figure gallery](../figures/osteosarc/2026-09-16_05-43-34-415330Z/README.md)
or the [original three-example, 12-page vector PDF](../figures/osteosarc/2026-09-16_05-23-24-826750Z/three-examples-12-pages.pdf).
To regenerate them from a repository checkout with the plotting extra installed:

```sh
python -m examples.osteosarc_assembly_figures --output-dir figures/osteosarc
```

This uses checksum-pinned, original, primary-only RNA fixtures and offline
Ensembl 87 models. No download, altered read, synthetic read or realignment is
involved. It writes twelve example directories, with their captions and a
manifest. Each RNA protein must pass independent translation, frame and interval
checks. Each Varcode prediction must match an independent translation of the
reference plus the edit. Whether the RNA protein equals that prediction is
checked separately, not assumed.

| Example | Assembly off → on | Mutation-overlapping 25-mers |
| --- | --- | --- |
| DIAPH1 | 40 → 49 aa | 16 → 25 |
| SLC25A12 | 44 → 49 aa | 20 → 25 |
| TECPR1 | 37 → 49 aa | 13 → 25 |

In these three examples no single read spans the whole reconstructed cDNA, so
assembly recovers context by joining overlapping reads. The shorter sequence
without assembly is correct, but misses those peptide windows. Fixture counts
are selections, not whole-sample abundance.

Other examples:

- **DYNC1H1:** the V314I and Q3267H loci, kept separate.
- **H1-2:** a 15-nt in-frame deletion; the pinned annotation names the gene
  HIST1H1C. Varcode's repeat-normalized protein deletion label differs from the
  source's, but its predicted sequence is checked independently.
- **MAP2:** a 22-nt frameshift deletion with one alternate RNA read, so no protein
  under the default coverage floor. Two reference predictions are shown, without
  inventing RNA support.
- **NAV2:** two transcripts give the 47-aa RNA sequence, which ends at W, before
  the reference-only tracks diverge into LR and VN. Eight annotation models are
  kept, seven with protein predictions. The cDNA crosses no junction, so the
  figure does not establish one isoform or justify adding either predicted
  C-terminal extension.
- **NTF3:** Varcode's nominated A>G predicts K56R (K69R in the other model), but
  the RNA encodes S. An adjacent change, linked to it in the source data, makes
  a compound AG>GT allele. The independent translation checks pass, and the single-edit
  comparison fails, as expected. This does not claim the second change is
  germline, or that assembly is needed: both modes recover S.
- **NR2F2:** native GRCh37 CeGaT RNA. Three selected original templates carry the
  mutation and a downstream 3-nt CIGAR deletion. Both assembly modes keep the
  deletion, and the single-edit Varcode prediction lacks it. The matched-sample
  audit finds the deletion in RNA but in none of the assessed tumor-DNA or
  blood-DNA templates, and the RNA templates share one alignment and
  fragment-endpoint family. Treat it as RNA-supported context not confirmed in
  DNA, not as a validated germline or somatic deletion.
- **PIP5K1A, T1 ONT versus T1 Illumina:** complete regional inputs rather than
  selected fixtures. ONT gives 46 aa in both modes; Illumina has one alternate
  read and no protein at the two-read floor. This is a difference in library
  support, **not a controlled read-length comparison**. See the
  [pinned original-read comparison and reproduction recipe](../tests/data/osteosarc/figure_comparisons/README.md).

### Structural variants

Literal indels, symbolic SV calls and supplied fusion transcripts are different
inputs. The small-variant pipeline rejects symbolic alleles
([#270](https://github.com/openvax/isovar/issues/270)) rather than inferring a
rearranged protein from `END` and a gene label. Use
[`isovar sv-rna`](sv-rna.md) for an SV call and [`isovar fusion`](fusion.md) for
a supplied fusion transcript.

The expanded gallery includes original fusion RNA and the
[nine-candidate RNA footprint audit](../tests/data/osteosarc/figure_comparisons/EXTENDED_RNA.md),
which queries SV loci directly, without needing a Varcode protein prediction.
DLG5 shows why DNA assembly, RNA splice footprints and translated consequences
must stay separate: its sequence-resolved DNA junction is supported, but no
mutant RNA or CDS is established. A novel splice path can be a candidate
consequence, but it is not assigned to the mutation without direct or safely
phased linkage. `python -m examples.osteosarc_context_figures` writes the full
dated gallery and a combined vector PDF.

The public data are CC0; see the
[source and selection documentation](../tests/data/osteosarc/README.md). Each
batch records the exact source BAM URLs, original alleles, fixture and reference
checksums, effective settings and independent checks, without read names.

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

- `collect_visualization_data` does not import Matplotlib. `creator_kwargs`
  passes `ProteinSequenceCreator` settings, and
  `include_reference_predictions=False` leaves out the Varcode tracks without
  changing the RNA analysis.
- The data include `reference_predictions`, `varcode_version` and, per
  transcript, `rna_supported`: whether the model contributes to a displayed
  protein, which does not identify a unique isoform.
- In each mode, `protein` is the top result and `proteins` every result allowed
  by the recorded cap. Every protein keeps its cDNA and frame contexts and one
  real reconstruction (its `witness`).
- `plot_variant_evidence(data)` returns a Matplotlib figure for your own
  styling. It does not select a GUI backend or change global plotting defaults.

To compare several sequencing products, collect each one separately, with the
same annotation and filters, and render them together:

```python
from isovar.protein_comparison import save_protein_comparison

products = []
for source_id, label, read_evidence in independently_collected_products:
    data = collect_visualization_data(
        variant, read_evidence, compare_assembly=True,
        creator_kwargs={"max_protein_sequences_per_variant": None},
    )
    products.append({"source": source_id, "label": label, "visualization": data})
save_protein_comparison(products, run / "all-products")
```

Products are never pooled into one count. Missing coverage, no alternate
support and failure to translate a protein are shown as different outcomes, and
reference sequence never fills missing RNA.
