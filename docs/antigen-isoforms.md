# CTA and other antigen isoforms

`isovar isoforms` reconstructs RNA at explicitly selected genes without a
nominated mutation. It keeps competing transcripts and sequence alleles,
reports ambiguous fragment assignments, and groups complete encoded protein
sequences. You supply the CTA target list: expression alone does not establish
cancer-testis specificity or expression in malignant cells.

## Reconstruct a panel

Use the FASTA used for alignment and a matching, versioned GTF. For new human
analyses, the preferred combined alignment reference is GDC **GRCh38.d1.vd1**;
see [reference annotation](reference-annotation.md#combined-human-and-viral-reference).
GENCODE's comprehensive primary-assembly annotation retains more isoforms than
the basic subset. The annotation version is an explicit input, separate from
the alignment reference.

```sh
isovar isoforms reconstruct \
    --bam tumor-rna.bam --sample-id tumor-1 \
    --gene PRAME --gene MAGEA4 --gene CTAG1B \
    --genome-fasta GRCh38.d1.vd1.fa \
    --gtf gencode.v50.primary_assembly.annotation.gtf.gz \
    --annotation-version 50 \
    --output-dir tumor-1-isoforms
```

These gene names illustrate an explicit panel, not an exhaustive CTA definition.
Include related loci when reads could map to both. BAM secondary placements,
mates and declared unreported alternatives remain part of the attribution.
`--alignment-reference custom --assembly GRCh37 --annotation-name Ensembl`
supports another assembly/reference; set its actual annotation version too.

The command writes:

- `isoforms.json`: sequences, coding interpretations, reference identities,
  per-base coverage, original read witnesses, and complete-input support.
- `transcripts.fa`: full annotated competitors and reconstructed paths reaching
  annotated transcript boundaries. A reference competitor can contain uncovered
  sequence; its origin and coverage are explicit in the JSON.
- `proteins.fa`: complete coding hypotheses with RNA coverage over every coding
  base. The JSON distinguishes a whole-CDS read witness from overlap assembly.

Observed partial sequences remain in `partial_sequences`; they are never filled
with reference bases to make full transcripts. `protein_groups` combines identical
complete amino-acid sequences while retaining all transcript identities and
counting each compatible fragment once. An encoded sequence is not proof of
translation or antigen presentation.

Each whole-read assembly seed retains its original placement choice through
both extension directions. Alternative alignments of the same sequenced
segment cannot extend one another into a longer RNA sequence. They remain
separate hypotheses; independently observed overlaps can still extend a path.
Discovery retains original parent identities even when sequences and placements
are identical, so sharing geometry cannot hide an independent overlap witness.

The default discovery budget is 10,000 target read segments and 256 assembled
paths. The later support pass checks the complete declared input, even when
discovery stops early. Increase `--max-discovery-reads` or `--max-paths` when the
result reports truncation. A truncated search cannot establish a dominant
isoform. Disk-backed indexing avoids retaining background alignments; the
in-memory report still grows with its retained target witnesses.
Annotation indexes are shared across samples in PyEnsembl's cache. Use
`--reference-cache` to choose another cache root; annotation file contents
select a separate child directory, preventing stale or same-basename indexes.
`support.fragments_without_full_transcript_assignment` also counts fragments
assigned only to partial sequence, which cannot be quantified as a full isoform.

## Shared read policy and library reuse

`ReadCollector.for_rna_paths()` supplies the same default policy for literal
variant paths, CTA/antigen isoforms and viral reconstruction: include secondary
placements, MAPQ-zero reads and missing qualities; exclude duplicate-marked
records and soft-clipped bases. Both discovery and full-input support use the
same collector, and `collection` records the effective settings/profile digest.

The reconstruction CLI accepts the common `--min-mapping-quality`,
`--use-duplicate-reads`, `--drop-secondary-alignments`, `--require-base-qualities`,
`--use-soft-clipped-bases`, `--infer-read-ends`, `--read-end-profile`,
`--trim-adapters` and `--trim-poly-a` options. Adapter profiles are explicit;
see [technical-sequence evidence](read-end-inference.md). Original read records
remain unchanged. Original mate/segment identities are retained; small-variant
mate-merging overrides are not controls for this path workflow and are rejected.

Use `--library-cache DIRECTORY` to reuse an original BAM index. The Python API
accepts `read_collector`, `library_index` (an open `RnaLibraryIndex`) and
`read_filter_identity` for a custom eligibility predicate. Build the index with
`prepare_rna_library(..., read_collector=collector, reference_name=...)` using
the same policy/reference. Reuse validates the original source, header,
reference context and eligibility; sequence policies are applied again to the
original reads. Ownership of a supplied index stays with the caller. Persistent
indexes require a local seekable BAM; uncached input also supports SAM/CRAM.

Cohort `viral` and `isoforms` tasks accept the same `read_collection` mapping,
including a profile from the manifest's verified file registry. Tasks sharing a
source, reference context and policy reuse an index within the cohort run.
Separate libraries keep separate identities. Put read options in
`read_collection`, not `parameters`; changed settings invalidate result reuse.

## Estimate abundance with Salmon

Reconstruction support counts are not transcript abundance. For short reads,
the optional second command runs installed Salmon against your **full competing
transcriptome**, replacing the selected reference entries with the reconstructed
competitors. Shared reads are not assigned to each transcript independently.

```sh
isovar isoforms quantify \
    --input tumor-1-isoforms/isoforms.json \
    --transcriptome gencode.v50.transcripts.fa.gz \
    --reads1 tumor_R1.fastq.gz --reads2 tumor_R2.fastq.gz \
    --threads 4 --output-dir tumor-1-abundance
```

Omit `--reads2` for unpaired short reads. This command is not a long-read
quantifier. Run samples/libraries separately; do not concatenate related signal
products and treat them as independent observations.

The output includes `isoforms.json`, the augmented FASTA, Salmon's index and
`quant/quant.sf`, logs, commands, version, and input checksums. Target transcript
IDs and sequences must agree with the supplied background. Existing output
directories are rejected. No program or reference is downloaded automatically.
The integration was exercised with Salmon 2.8.0.

`quantification.genes` reports within-gene fractions and the top estimated
isoform when it has distinguishing fragment evidence and discovery completed.
Indistinguishable transcript groups remain unresolved even if Salmon assigns
them different values. The top estimate does not assert statistically significant
dominance; no inferential confidence interval is computed by this wrapper.
`quantification.proteins` sums transcript TPM for identical complete proteins,
separating RNA-covered coding sequence from uncovered reference hypotheses.
Alternative unresolved coding interpretations are not each given the full
transcript abundance.

The denominator is the supplied competitive transcriptome. A regional panel
alone cannot establish whole-transcriptome TPM. Keep background paralogs and
other transcripts in the FASTA, and use reads from the declared sample/library.
The reconstruction and quantification input identities are recorded separately.

## Python

```python
import pysam
from isovar import (
    gdc_genome, isoform_targets_from_genome, validate_alignment_reference,
    reconstruct_antigen_isoforms, quantify_antigen_isoforms,
)

with gdc_genome(
    genome_fasta="GRCh38.d1.vd1.fa",
    gtf="gencode.v50.primary_assembly.annotation.gtf.gz",
    annotation_version=50,
    cache_directory="annotation-cache",
) as genome, pysam.AlignmentFile("tumor.bam") as bam:
    genome.index()
    targets, unavailable = isoform_targets_from_genome(
        genome, ["PRAME", "MAGEA4"], contig_names=bam.references)
    reference = validate_alignment_reference(
        bam, genome, profile="GRCh38.d1.vd1",
        contigs=[t.reference.contig for t in targets])
    result = reconstruct_antigen_isoforms(
        bam, targets=targets, sample_id="tumor-1", source="RNA-library-1",
        reference_identity=reference)
    result["target_selection"] = dict(
        genes=["PRAME", "MAGEA4"], unavailable_transcripts=unavailable)

quantify_antigen_isoforms(
    result, transcriptome="gencode.v50.transcripts.fa.gz",
    reads1="tumor_R1.fastq.gz", reads2="tumor_R2.fastq.gz",
    output_dir="tumor-1-abundance")
```

`IsoformTarget` wraps the existing `FusionReference` sequence/exon/CDS model
with gene identity. Callers with explicit local/viral models can construct it
directly. CDS end includes the stop codon; exon intervals are zero-based and
half-open. Coding hypotheses currently use the standard genetic code and an
annotated start/frame. Unsupported reference CDS models are retained without
inventing an ORF; ambiguous reference bases make that transcript unavailable.

## Current boundaries

This first part of [#465](https://github.com/openvax/isovar/issues/465) uses exact
coordinate-anchored overlaps. It preserves observed sequence differences; it
does not decide whether each difference is an allele or a sequencing error.
Noisy long reads can therefore remain unassigned or yield additional hypotheses.
Overlap assembly is not end-to-end phasing, and reaching annotated ends does not
prove a physical RNA molecule's transcription start or polyadenylation site.

[`isovar viral`](viral-sequences.md) recovers observed viral nucleotide paths
from GDC alignments without a GTF. GENCODE does not annotate those viral
transcripts; annotated transcript/protein interpretation needs versioned,
sequence-matched models. A composite GTF can use `--annotation-name` and its own
`--annotation-version`; overlapping CDS, programmed frameshifts, polycistronic
translation and mature polyprotein processing need further explicit support.
Their general handling, noisy-read inference, and a validated real multi-sample
CTA/viral corpus remain tracked in #465.

Methods: [GENCODE formats](https://www.gencodegenes.org/pages/data_format.html),
[Salmon's primary paper](https://www.nature.com/articles/nmeth.4197), and
[Salmon CLI documentation](https://combine-lab.github.io/salmon/reference/cli/).
