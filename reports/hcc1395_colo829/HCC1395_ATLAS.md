# HCC1395 whole-run original-RNA atlas

This report joins the [SEQC2 v1.2.1](https://pmc.ncbi.nlm.nih.gov/articles/PMC8532138/)
HCC1395/HCC1395BL matched-normal SNV/indel benchmark to original ONT RNA
from [PRJNA934933](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA934933).
The [Talsania et al. 2022](https://link.springer.com/article/10.1186/s13059-022-02816-6)
consensus SV set is screened separately. The RNA libraries and DNA benchmark
were made in different studies, so their cell passages cannot be equated.
An RNA reference observation does not negate a SEQC2 DNA somatic call.

## Inputs and method

- The source universe contains 39,447 SEQC2 SNVs and 1,625 indels: 38,797
  `PASS;HighConf` and 2,275 `PASS;MedConf`. Both confidence classes are
  enumerated; they are not pooled into a single high-confidence label.
  HCC1395BL is the normal DNA baseline. Every one of the 41,072 REF alleles
  matches the NCBI GRCh38 no-alt FASTA at its one-based VCF position.
  The source VCF SHA-256 values are in [README.md](README.md#input-sha-256).
- The Talsania Table S4 source contains 1,788 consensus SV records: 717 DEL,
  551 INS, 230 DUP, 146 TRA, 133 INV, and 11 BND. Table S9 names 24
  RNA-seq fusion associations; 20 IDs match a Table S4 record. These are
  published external annotations, separate from this original ONT run.
- ONT direct RNA [SRR23456709](https://www.ncbi.nlm.nih.gov/sra/SRR23456709)
  has 1,049,452 reads and 1,068,238,021 bases according to SRA Toolkit 3.4.1
  `sra-stat --quick`. The validated 1,006,368,588-byte SRA archive has
  SHA-256 `d69db8e48bbb20665c3e03a1fa35742be818a0127994e015cb3a560afae39ac8`.
  The [ENA run report](https://www.ebi.ac.uk/ena/portal/api/filereport?accession=SRR23456709&result=read_run&fields=run_accession,sample_accession,study_accession,library_strategy,library_source,library_layout,read_count,base_count,fastq_bytes,fastq_md5&format=tsv)
  assigns BioSample SAMN33292477 and FASTQ MD5
  `1f5d2a60c3158ae12e573a316bd846bb`.
- All original reads were streamed through minimap2 2.26
  (`-ax splice -uf -k14`) against the NCBI GRCh38 no-alt analysis set,
  uncompressed FASTA SHA-256
  `9cce8b926416dd96b152deea85188495b75f7ac8d634cc723a017067be8702b7`.
  An indexed BAM retains reads overlapping any SEQC2 allele or either 500-base
  SV breakend window; this is output selection **after the complete run was
  aligned**. Its SHA-256 is
  `9634b78414b9ddbbb9e098a01d899028dd0b2d703c3c9d0a5c24ce6a95cc9db9`.
  The selected BAM contains 90,603 alignment records, including 80,658
  primary and 685 supplementary records. It passes `samtools quickcheck`.
- Isovar 1.59.7 with Ensembl 87 counts allele-read objects after its default
  read handling; these single-end runs have one original read per fragment.
  Translation uses the Isovar default local context and coverage rules.
  Protein rows retain every ranked context, not a single forced isoform.

## Direct RNA small-variant findings

All 41,072 SEQC2 small-variant records have an RNA and translation status.
The 343 ALT-observed calls comprise 333 SNVs and 10 indels; 339 have the
source `PASS;HighConf` label and four have `PASS;MedConf`. They account for
2,681 ALT allele-read observations, with original read names in the witness
table. Local translation yielded 119 ranked protein contexts across 22
calls, including one indel. The remaining 321 ALT-observed calls yielded no
local context under these settings.

| RNA status | Calls |
|---|---:|
| `alternate_observed` | 343 |
| `reference_only` | 239 |
| `reference_and_other_only` | 53 |
| `other_allele_only` | 19 |
| `alignment_without_allele_call` | 7,300 |
| `no_alignment_overlap` | 33,118 |

`alignment_without_allele_call` includes alignments whose reference span
overlaps a site without a callable base at that site, such as a skipped
intron. It is not allele coverage. A no-overlap status applies only to the
specified run, alignment, and extraction method; it does not demonstrate
absence of expression.

Two independently reported mutant-peptide loci gain direct RNA local
contexts in the full run:

| SEQC2 call | DNA tumor/normal VAF | RNA REF/ALT/other | Top local mutant context | Zero-based mutant interval | Supporting reads |
|---|---|---:|---|---:|---:|
| OSTC F9L, `COSM32736` | 0.490 / 0.001 | 75 / 60 / 10 | `METLYRVPLLVLECPNLKLKKPPWL` | `[8,9)` | 5 |
| PRDX5 F157L, `COSM32759` | 0.442 / 0.002 | 362 / 193 / 155 | `VRLLADPTGALGKETDLLLD` | `[10,11)` | 6 |

The DNA VAF values are the SEQC2 VCF `TVAF`/`NVAF` annotations, not measured
RNA fractions. [Independent proteogenomics](https://pmc.ncbi.nlm.nih.gov/articles/PMC11419083/)
reports mutant residues in OSTC and PRDX5. This supports the residue
comparison, not these full local amino-acid contexts, full-length isoforms,
protein abundance, or cell-passage identity. The 100,000-read subset in
[README.md](README.md) had ALT reads at both loci but lacked a default direct
RNA translation; the complete run resolves local contexts with five and six
supporting reads. A SEQC2 indel at `chr8:123220221:A>AGGGCTAC` has three ALT
reads and a top local context `HPEVKCVGSVARSPDCLGDCILRRI` with two
supporting reads.

## Direct RNA structural-variant screen

Each of the 1,788 Talsania SV rows has observed original-read counts in both
500-base breakend windows, plus a same-read split-alignment candidate count
when the two windows are distinct. Overlapping windows are not compared for
split candidates. The outcomes are:

| RNA status | SV records |
|---|---:|
| `no_breakend_window_alignment` | 1,233 |
| `one_end_covered_no_split_candidate` | 149 |
| `both_ends_covered_no_split_candidate` | 400 |
| `split_alignment_candidate_unresolved` | 6 |

The six candidate records have 17 record-level split-read observations in
total; one read may contribute to multiple nominated DNA records. The
published POU2F1-NME7 association has one observed same-read split candidate
(`SRR23456709.320596`), but this screen has not checked junction bases,
orientation, splice interpretation, reading frame, or a translated fusion
ORF. Its `path_status` and `protein_status` therefore remain unresolved.

The separate [junction audit](sv_junction_audit.tsv) replays three oriented
events against a source-derived original-read BAM, retaining low-MAPQ pieces
and soft-clipped sequence. Isovar finds 7 direct junction fragments for the
chr9 duplication, 4 for the chr13 deletion, and 2 for the chr21 duplication;
the MAPQ 20 sensitivity subset retains 7, 3, and 2. Isovar extends the
junctions through overlapping original read sequences; one read need not
cover the resulting context. The chr9 search reaches its ten-path cap.
These are RNA junction observations and bounded
assembly hypotheses, not unique full-length isoforms or SV-derived protein
calls. The three chr1 candidate records reuse one complex read and remain
outside this oriented replay. The broad screen's unresolved statuses describe
its original 500-base proximity test; the junction audit is a deeper,
selected-event result.

## Outputs and replay

- `hcc1395_direct_rna_atlas_calls.tsv.gz` assigns every SEQC2 small variant
  an allele status, read/fragment counts, genomic transcript overlaps, and
  translation status. Genomic overlaps are reference annotations, not
  observed RNA-compatible isoforms.
- `hcc1395_direct_rna_atlas_proteins.tsv.gz` contains all ranked, locally
  translated contexts with mutation interval, compatible Ensembl transcript
  IDs, and original read/fragment support. They are not resolved full-length
  mutant isoforms. A noisy ONT ALT observation alone does not establish one.
- `hcc1395_direct_rna_atlas_alt_witnesses.tsv.gz` retains original SRA read
  names and source-read counts per ALT observation.
- `hcc1395_direct_rna_sv_screen.tsv` retains all source SV records, both
  breakend counts, candidate read IDs, and explicit unresolved path/protein
  statuses. Published fusion names are copied from Table S9 for comparison.
- `tests/data/hcc1395_whole_run/` has an indexed, 1.4 MB original-read BAM
  around informative calls and SV ends for offline regression tests. The
  full selected-site BAM is 149 MB and stays out of Git and CI.

Rebuild the full alignment after downloading the source archive and FASTA
named above:

```sh
python -m reports.hcc1395_colo829.build_hcc1395_regions \
  --output hcc1395-regions.bed
bash reports/hcc1395_colo829/align_full_hcc1395_rna.sh \
  GRCh38-no-alt.fa hcc1395-regions.bed SRR23456709.sra \
  /path/to/sratoolkit/bin/fastq-dump atlas-work direct_rna
python -m reports.hcc1395_colo829.build_hcc1395_atlas \
  --bam atlas-work/SRR23456709.all_sites.grch38.bam
python -m reports.hcc1395_colo829.build_hcc1395_sv_screen \
  --bam atlas-work/SRR23456709.all_sites.grch38.bam
```

The archive and full selected-site BAM are not required for the bounded
offline tests. The atlas tables are whole-run outputs and should not be
regenerated from the bounded fixture. The DNA matched-normal calls,
direct RNA read witnesses, local translations, published fusion associations,
and unresolved SV paths remain separate evidence types.
