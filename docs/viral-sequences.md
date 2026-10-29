# Recover viral sequences from combined-reference RNA alignments

When RNA has been aligned to GDC **GRCh38.d1.vd1**, Isovar can recover observed
viral nucleotide paths directly. No VCF, gene list or human GTF is required:

```sh
isovar viral \
    --bam tumor-rna.bam \
    --genome-fasta GRCh38.d1.vd1.fa \
    --sample-id tumor-1 \
    --output-dir tumor-1-viral
```

This examines the 200 viral contigs explicitly listed by GDC, including `chrEBV`,
`HPV16`, `HPV18`, `HBV`, `HCV-1`, `KSHV` and `MCV`. It verifies the full BAM/FASTA
sequence dictionary and the selected viral DNA against GDC's published MD5s.
The panel's source URL, checksum and historical reference IDs are packaged and
recorded in each result. The reference FASTA remains a supplied local input.

Use `--contig HPV16 --contig chrEBV` to examine a subset. For another alignment
reference, use `--alignment-reference custom` and explicitly supply its viral
contigs with `--contig`. Custom references are checked against the BAM's lengths
and available M5 sequence checksums; missing M5 leaves sequence identity
explicitly unverified. A human-only BAM cannot supply viral alignments: retain
the original sequences and align competitively to the combined reference first.

The output directory contains:

- `sequences.fa`: observed nucleotide paths, including partial reconstructions.
- `proteins.fa`: predicted complete and partial coding products anchored to
  available ORFdata references. Headers retain RNA sequence IDs, coding status,
  orientation and reference ORF anchors. This file is empty when no coding
  products can be annotated.
- `viral.json`: schema `isovar.viral_sequences.v1`, with sequences, mapping runs,
  per-base coverage, contig counts, original fragment/segment identities,
  alignment metadata, barcode/UMI and lineage evidence, discovery limits and
  source/reference checksums, coding predictions and coding-reference coverage.

Sequence and splice differences remain separate. Uncovered bases are never
inserted from the viral genome, and an entire viral genome is not treated as a
transcript. Paths containing an observed host/viral split retain their host
sequence too, labelled by `other_contigs` and the original placements. Such a
path is evidence of a chimeric alignment, not by itself proof of integration.

The canonical sequence orientation supports overlap assembly but does not infer
RNA strand from an unstranded library. `rna_strand` and
`transcript_completeness` remain `unknown`. A complete read witness and complete
coverage of a reconstructed path are reported independently of transcript ends.
Noisy reads can remain unassigned or create extra exact-sequence hypotheses.

## Counts and resource limits

`contigs` reports aligned fragments/segments, unambiguous viral-contig fragments,
sequence-assigned fragments and unassigned fragments. The support pass considers
the complete declared input, retaining mates and secondary/supplementary records.
RG/QNAME identify a fragment; mate bits distinguish segments. Identical original
SAM payloads count once. A fragment with placements on two viruses contributes
to both contig totals, so summing contig totals can exceed the number of original
fragments. The report's global `support.target_fragments` counts each once.
An unambiguous viral-contig assignment can still leave several sequence
placements on that same virus; the sequence-level `unique_fragments` counts
retain that distinction.

Shared paths remain ambiguity groups. Missing mapped mates, missing/malformed SA
declarations, unreported NH placements and ambiguous/unavailable MAPQ prevent
unique attribution. The report keeps their original declarations and explicit
uncertainty reasons. Supplied barcode and UMI labels do not establish independent
molecules or expression in malignant cells. These fragment counts are neither
TPM nor viral genome copies.

By default, discovery uses up to 10,000 observed segments and 256 paths **per
viral contig**, with a 30-base minimum exact overlap. Each contig has its own
budget, and discovery retains only one contig's observations at a time. A disk
index keeps background records and the original mates/placements available.
Use `--scratch-dir` to choose temporary storage and `--reference-cache` to choose
PyEnsembl's genomic-FASTA cache. Support evidence and output memory still grow
with the retained target witnesses and reconstructed paths.

Increase `--max-discovery-reads` or `--max-paths` when a contig reports
`discovery_truncated`. Later reads still count toward support even after the
discovery budget is exhausted. Contigs with no eligible alignments report
`no_eligible_alignments`, not viral absence; excluded-record counts and
unassemblable aligned reads remain explicit. Reference strain differences,
discarded/unmapped reads and unavailable placements cannot be recovered from a
filtered BAM by adding this command.

## Shared policy and indexes

Viral reconstruction uses the same [read-policy and index controls](antigen-isoforms.md#shared-read-policy-and-library-reuse) as antigen isoforms and literal RNA paths. The API accepts `read_collector`, `library_index` and `read_filter_identity`; the CLI exposes common read/adapter flags and `--library-cache`. Cohort `read_collection` applies to both discovery and support. Effective settings are recorded in `collection`.

## Coding products and ORFdata

Isovar 1.57 requires ORFdata 0.29 or newer. ORFdata owns the small, offline coding
catalogue and its provenance; Isovar does not duplicate CDS annotations or add
viral transcripts to human Ensembl releases. The initial catalogue has the nine
CDSs annotated in RefSeq **NC_001526.4**, including its two joined products. It
recognizes both that native genome and the exactly matched, rotated GDC `HPV16`
sequence by their full DNA checksum, independently of contig names.

This is a source CDS subset, not an exhaustive HPV16 transcript or ORF catalogue.
Other GDC viruses still produce RNA sequences and counts, with explicit
`unavailable_reference` coding status. A differently sequenced HPV16 reference
also remains unavailable rather than inheriting annotations by name. Broader
coverage and special coding semantics are specified in
[ORFdata's viral roadmap](https://github.com/openvax/orfdata/issues/91), with
separate issues for [other HPV types](https://github.com/openvax/orfdata/issues/93),
[EBV/KSHV](https://github.com/openvax/orfdata/issues/94),
[CMV](https://github.com/openvax/orfdata/issues/95),
[hepatitis viruses](https://github.com/openvax/orfdata/issues/96),
[HTLV-1/HIV](https://github.com/openvax/orfdata/issues/97), and
[polyomaviruses](https://github.com/openvax/orfdata/issues/98).

For each observed path, Isovar uses a placed reference codon to anchor a reading
frame and translates only observed RNA bases to the first observed stop or the
last complete codon. Both orientations are considered without asserting the
library's RNA strand. Each `coding` entry records query coordinates in the
original RNA, coding sequence, amino acids, observed start/stop status and
reference ORF associations. Query coordinates are 0-based, half-open; reverse
predictions still use the original RNA's coordinates. Shared predictions retain all their reference
anchors. Altered splice/indel paths and premature or displaced stops remain
explicit. A reference anchor does not assert that the sample encodes that
reference's original protein; `protein_matches_reference` compares complete
predictions and is unknown for partial products.

`complete_CDS` means the source start and an in-frame stop are observed in this
prediction. It does not establish transcript boundaries, initiation, protein
expression or peptide presentation. `partial_CDS` retains missing coding ends;
unobserved bases are never supplied from the reference. Coding support reports
observed base coverage and whole-CDS segment witnesses independently of RNA
fragment counts. Neither predicted proteins nor source annotations supply
experimental translation evidence.

The top-level `coding_reference` records the ORFdata version, per-contig DNA
checksums and coverage. Available references retain accession/source checksums,
catalogue checksum, source-coordinate rotation and product inventory, allowing
the annotation set to be identified even when the installed catalogue changes.

## Python API

```python
import pysam
from isovar import reconstruct_viral_sequences, write_viral_fasta, write_viral_proteins

with pysam.AlignmentFile("tumor-rna.bam") as bam:
    result = reconstruct_viral_sequences(
        bam, genome_fasta="GRCh38.d1.vd1.fa",
        sample_id="tumor-1", source="RNA-library-1")
write_viral_fasta(result, "viral-sequences.fa")
write_viral_proteins(result, "viral-proteins.fa")
```

PyEnsembl's ordinary `Genome` provides genomic DNA and its FASTA cache without
requiring a GTF. ORFdata provides sequence-qualified `SourceRNA` and
`TranscriptORF` models independently of human annotations. Translation exceptions,
programmed recoding and polyprotein processing require the corresponding
upstream models; this release does not infer them from an arbitrary RNA path.

## Recovering reads that lack viral alignments

This command consumes alignments. For paired short RNA reads, the recovery
experiment in [#514](https://github.com/openvax/isovar/issues/514) favored direct
splice-aware minimap2 alignment to a small viral panel with
`-ax splice:sr -K 1m -t 4 --sam-hit-only`. Small query batches retained the tested
alignments while using less memory than the default batch; adding a separate
k-mer filter did not improve the measured runtime/sensitivity tradeoff. No
screen is guaranteed to recover every viral read.

Panel-only hits are candidates. Retain original mate sequences, qualities,
RG/QNAME identities and reported alternatives, then align candidates
competitively against human/decoy/viral DNA before interpreting viral support.
The full combined reference is downloadable from the GDC sources below. For
memory-limited verification, minimap2 supports `-I 100M --split-prefix` with a
temporary prefix to merge split-reference results before computing MAPQ.
`-K 1m` bounds query batching; it does not bound the full reference index. These
are external alignment steps, not additional `isovar viral` options. See the
[minimap2 manual](https://github.com/lh3/minimap2/blob/v2.31/minimap2.1).

Sources: [GDC reference files](https://gdc.cancer.gov/about-data/gdc-data-processing/gdc-reference-files)
and its [viral sequence list](https://gdc.cancer.gov/system/files/public/file/GRCh83.d1.vd1_virus_decoy.txt).
