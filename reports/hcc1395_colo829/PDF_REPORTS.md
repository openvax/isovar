# HCC1395 and COLO829 somatic RNA evidence PDFs

The [HCC1395 PDF](../../output/pdf/HCC1395_somatic_RNA_evidence.pdf) and
[COLO829 PDF](../../output/pdf/COLO829_somatic_RNA_evidence.pdf) summarize the
whole-run atlas tables for the selected public matched-normal DNA sources.
They include every small-variant record in those sources and every selected
SV source record. They are not an exhaustive list of mutations in every
culture of either cell line. The per-call results remain in the compressed
atlas tables; the PDFs give counts, representative evidence, methods, and
limits.

| Measure | HCC1395 | COLO829 |
|---|---:|---:|
| Ingested small variants | 41,072 | 46,064 |
| Ingested SV source records | 1,788 | 188 |
| SV source records involving an Ensembl 87 exon | 354 | 97 |
| Any Ensembl 87 exon involved | 1,728 | 1,395 |
| Exonic with original RNA ALT | 292 | 502 |
| Exonic with a local mutant protein context | 22 | 124 |
| Ensembl 87 coding-sequence effect | 491 | 373 |
| Pinned orfdata ncORF coding-block overlap | 17 | 5 |
| Either coding annotation | 505 | 378 |
| Original RNA ALT-observed calls | 343 | 1,933 |
| ALT-observed calls with a local mutant protein context | 22 | 124 |
| Ranked local mutant protein contexts | 119 | 12,218 |
| DNA-predicted top-effect frameshifts | 9 | 5 |
| Resolved SV-derived proteins | 0 | 0 |

The first figure in each report is a **nested source-record funnel**: all
ingested DNA records (small variants plus nominated SV source records), the
small-variant subset, any annotated exon involved, exonic with an original RNA
ALT allele, and exonic with a translated local protein context. An exon hit means
some changed reference base meets any Ensembl 87 exon, whether coding, UTR or
noncoding. A deletion is queried across all removed reference bases, so an
exon between intronic endpoints counts. Pure insertions have no removed base;
either flank is a *candidate* exon involvement. This follows the
[VCF padding-base convention](https://github.com/samtools/hts-specs/blob/master/VCFv4.3.tex)
and PyEnsembl's 1-based inclusive exon intervals. SV source records leave the
funnel at its second stage; breakpoint RNA evidence is screened separately.

Figures show ALT-read support separately for SNVs, insertions and deletions,
first with zero-support calls and then conditional on at least one ALT read.
The next page repeats those distributions for the **exon-involved subset** of
the same small-variant records. A separate figure groups exon-involved
nominated DEL, DUP, INV, INS, TRA and BND records by RNA split-alignment
candidate count, again with and without zeros. A structural variant's full
nominated intrachromosomal span is used for DEL/DUP/INV; INS uses flanks and
TRA/BND use nominated breakpoints. Four COLO829 breakend records cannot be
ruled out as exonic because one partner coordinate is unresolved, and are
marked `Unknown` in the
[SV exon annotation table](sv_exon_annotations.tsv.gz). A positive split
candidate is a screen, not a validated event-linked junction or protein.
The SV annotation [metadata](sv_exon_annotations.tsv.metadata.json) pins the
two source-table SHA-256 digests and reference release. The figures do not
combine ALT-spanning reads with split-alignment candidates.

The changed-amino-acid figure counts positions in the top-ranked **local**
protein window of each ALT-observed call, retaining the no-product and
zero-length deletion-junction categories. Its second panel excludes the
one-amino-acid category to show the tail. It does not estimate an unobserved
frameshift tail or full-length isoform.

Example tables rank original RNA ALT support, highly supported non-SNV local
products and the least reference-self-like local windows. The last ranking is
an exact screen: for each top local product, count distinct
mutation-overlapping 9-mers absent from the pinned Ensembl 87 GRCh38 protein
FASTA. Rank by absent fraction, then absent count, then ALT reads. It is not a
patient-specific self inventory, peptide presentation or immunogenicity score.
The FASTA gzip SHA-256 is
`d873c199d5dbf20869db3f2cbd6ac726259215d91c2782f9f1ddcad4bb117cf9`.
The per-call scores are in
[the compact score table](reference_self_9mer_scores.tsv.gz); the source
[Ensembl 87 protein FASTA](https://ftp.ensembl.org/pub/release-87/fasta/homo_sapiens/pep/Homo_sapiens.GRCh38.pep.all.fa.gz)
is not bundled. DNA-predicted frameshifts, reference ncORF coding-block
candidates, and unresolved SV screens are kept separate from RNA-resolved
protein products.

The [per-call reference annotation table](reference_coding_orfdata_annotations.tsv.gz)
contains all 87,136 small-variant records. `exon_overlap=True` applies the
funnel's any-exon rule. `reference_coding_effect=True`
means at least one Varcode effect modifies an Ensembl 87 protein-coding
sequence. `orfdata_coding_block_overlap=True` means the reference event span
overlaps a coding block in the pinned orfdata 0.34.1 human GRCh38 catalog
(7,264 models, release 2026-10-11.1, SHA-256
`06bc0bf3f81c25a7356f71a38581ff5704c05dcb75752454bc1ca75075ab8c09`).
For an insertion, the VCF anchor is only a *candidate* overlap; neither it
nor another ncORF overlap proves a mutant ncORF translation. The three
HCC1395 calls shared by the two annotation categories are counted once in
the union. The COLO829 categories have no overlap. Five HCC1395 and two
COLO829 ncORF-overlap calls have an observed original RNA ALT allele; no
mutated ncORF protein was reconstructed from that observation.

The source RNA runs are HCC1395 ONT direct RNA SRR23456709 and COLO829
Illumina paired RNA SRR8615617. They are separate assays and cultures from
the matched-normal DNA call sets. An ALT allele read, local read-compatible
protein context, polished transcript sequence, and full-length isoform read
support are distinct evidence. Neither atlas resolves SV junction sequences
or calls an SV-derived protein. Viral RNA and targeted antigen-gene isoforms
were not assessed in these runs, so there is no viral distribution or
tumor-specific CTA protein count in these two PDFs. See
[HCC1395_ATLAS.md](HCC1395_ATLAS.md) and
[COLO829_ATLAS.md](COLO829_ATLAS.md) for the input hashes, alignment and
translation commands, source-specific counts, and per-call tables.

## Rebuild

The committed annotation table is generated by
`python -m reports.hcc1395_colo829.build_reference_annotations` with Varcode,
PyEnsembl release 87 data, and orfdata 0.34.1. The script checks the bundled
catalog SHA-256 and writes a metadata JSON sidecar. The optimized SNV path
uses the atlas's Ensembl 87 genomic transcript-overlap flag to skip Varcode
queries at sites without any overlapping transcript; every indel is queried.
The SV exon table is generated by
`python -m reports.hcc1395_colo829.build_sv_exon_annotations` with the same
Ensembl 87 annotation. This step needs the installed reference database;
rebuilding the PDFs from the committed tables does not.

Run `python -m reports.hcc1395_colo829.build_evidence_pdfs` with ReportLab
installed to regenerate both PDFs in `output/pdf`. The builder reads only
committed atlas and annotation tables; it does not download genomes, BAMs, or
call sets. It produces deterministic PDF bytes. The two offline regression
tests check the annotation tables against every small-variant and SV source
record, as well as selected known ncORF hits.

To regenerate the reference-self score table, supply the pinned Ensembl 87
protein FASTA to `python -m reports.hcc1395_colo829.build_reference_self_scores
--fasta PATH`. The script checks its SHA-256. The PDFs then rebuild without
the reference FASTA.

## Run on another variant/RNA set

Install `isovar[report]`, then run the installed command on a literal
small-variant VCF and an indexed tumor RNA BAM/CRAM:

```sh
isovar evidence-report --vcf somatic.vcf.gz --bam tumor-rna.bam \
    --sample-id tumor-1 --output tumor-1-evidence.pdf
```

This writes `tumor-1-evidence.pdf` and `tumor-1-evidence.tsv` with every loaded
variant, including zero-support records. `--nc-orf-catalog PATH` adds candidate
reference ncORF coding-block overlaps, after a species/assembly check.
`--self-fasta PATH` enables the exact reference-protein 9-mer screen. Without
those inputs, their categories are marked unassessed. The PDF has exon-involved
ALT-read plots. The VCF/BAM command does not accept symbolic SVs or infer
matched-normal somatic status; the caller must supply the appropriate DNA call
set. It uses the same high-depth direct-read context path as the COLO829 atlas
when more than 750 ALT reads occur. To compose separate SV, fusion, viral and
targeted antigen-isoform evidence from their versioned exports, see the
[antigen-source report command](../../docs/antigen-evidence-report.md).
