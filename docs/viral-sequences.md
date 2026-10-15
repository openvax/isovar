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
- `viral.json`: schema `isovar.viral_sequences.v1`, with sequences, mapping runs,
  per-base coverage, contig counts, original fragment/segment identities,
  alignment metadata, barcode/UMI and lineage evidence, discovery limits and
  source/reference checksums.

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

## Python and protein interpretation

```python
import pysam
from isovar import reconstruct_viral_sequences, write_viral_fasta

with pysam.AlignmentFile("tumor-rna.bam") as bam:
    result = reconstruct_viral_sequences(
        bam, genome_fasta="GRCh38.d1.vd1.fa",
        sample_id="tumor-1", source="RNA-library-1")
write_viral_fasta(result, "viral-sequences.fa")
```

PyEnsembl's ordinary `Genome` provides genomic DNA and its FASTA cache without
requiring a GTF. This command recovers RNA sequence; viral gene names, annotated
transcript models and encoded proteins need sequence-matched viral annotations.
Those can use the existing [isoform API](antigen-isoforms.md) when they fit its
supported coding model. General frameshifts, overlapping/polycistronic ORFs,
translation exceptions and polyprotein processing remain tracked in
[#465](https://github.com/openvax/isovar/issues/465).

Sources: [GDC reference files](https://gdc.cancer.gov/about-data/gdc-data-processing/gdc-reference-files)
and its [viral sequence list](https://gdc.cancer.gov/system/files/public/file/GRCh83.d1.vd1_virus_decoy.txt).
