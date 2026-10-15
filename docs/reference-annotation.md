# Choosing the reference annotation

Isovar takes transcript models, reading frames and gene names from a
[PyEnsembl](https://github.com/openvax/pyensembl) annotation. Install the
Ensembl release that matches your alignments before the first run, for example:

```sh
pyensembl install --release 75 --species human
```

Align the RNA to the same assembly as the annotation.

## On the command line

`--genome` names an assembly such as `GRCh38` or `hg19`, and the most recent
Ensembl release installed for it is used. To choose the release, add it after
the assembly: `--genome GRCh38:93` or `GRCh38.93` (as in Ensembl's file names)
uses Ensembl 93. A release that does not provide the assembly, such as
`GRCh38:75`, is an error.

## From Python

Pass the same name (`genome="GRCh38:93"`) or a PyEnsembl genome when loading
the variants:

```python
import pyensembl
import varcode
from isovar import run_isovar

variants = varcode.load_vcf(
    "cancer-mutations.vcf", genome=pyensembl.EnsemblRelease(93))
isovar_results = run_isovar(variants=variants, alignment_file="tumor-rna.bam")
```

A `pyensembl.Genome` built from your own GTF and transcript FASTA files works
the same way, and gene and transcript names then come from that annotation.

## Combined human and viral reference

For **new human alignments**, Isovar recommends GDC's **GRCh38.d1.vd1** combined
reference. It contains the GRCh38 no-alt analysis set, decoys and a viral panel.
It is the default checked profile of `isovar isoforms reconstruct`. Existing
`run`/VCF genome selection remains explicit: a new preferred reference cannot
change the coordinates or identity of an existing BAM.

Download the FASTA from [GDC reference files](https://gdc.cancer.gov/about-data/gdc-data-processing/gdc-reference-files).
The published archive MD5 is `3ffbcfe2d05d43206f57f81ebb251dc9`. Its viral panel
is an identified snapshot, not a catalogue of every current viral strain.

Keep three identities separate:

| Identity | Example | Purpose |
|---|---|---|
| Human assembly | GRCh38 | Human coordinates |
| Alignment reference | GRCh38.d1.vd1 | Exact human, decoy and viral sequences |
| Transcript annotation | GENCODE 50 | Exons, transcript IDs and coding frames |

For new CTA analyses, use the comprehensive GENCODE primary-assembly GTF for
the selected release. GDC's published STAR indexes instead use their stated
older annotation (v36 or v22); reproducing those alignments requires those exact
inputs. Using a newer annotation for interpretation does not retroactively
change the alignment's annotation. Missing patch/alternate-contig transcript
models are reported rather than silently mapped onto a primary chromosome.

PyEnsembl's existing custom `Genome` supports this arrangement. The convenience
factory adds no new reference type or annotation parser:

```python
from isovar import gdc_genome

genome = gdc_genome(
    genome_fasta="GRCh38.d1.vd1.fa",
    gtf="gencode.v50.primary_assembly.annotation.gtf.gz",
    annotation_version=50,
    transcript_fasta="gencode.v50.transcripts.fa.gz",
    protein_fasta="gencode.v50.pc_translations.fa.gz",
    cache_directory="reference-cache",
)
genome.index()
# This is a normal pyensembl.Genome; pass it to varcode.load_vcf(genome=...).
```

The transcript/protein FASTAs support the existing variant workflow. The new
[antigen isoform workflow](antigen-isoforms.md) derives its reference transcripts
directly from GTF exons and genomic DNA, so those FASTAs are optional there.
Always pair annotation and transcript/protein FASTAs from the same release.

`validate_alignment_reference` checks GDC's packaged dictionary of 2,779 contigs,
BAM/FASTA lengths, any BAM M5 checksums, and the DNA MD5 of requested target
contigs. The dictionary's source archive and digest are recorded with the
package. It does not claim to have sequence-verified unrequested contigs. For a
custom reference without BAM M5, sequence identity remains explicitly unverified.

Viral contigs retain their GDC names (for example `chrEBV` and `HPV16`); do not
rename a different viral accession to one of these names without checking the
sequence. Viral transcript/CDS annotation must be supplied separately. A
human-only BAM cannot establish viral absence, and adding an annotation cannot
recover reads discarded before that BAM was written.
