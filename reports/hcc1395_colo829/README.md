# HCC1395 and COLO829: source-call and RNA-isoform audit

Date: 2026-10-20, with 2026-10-21 whole-run extensions. This
inventory enumerates every record in the selected public DNA call sets and
every row in TAGET's COLO829 polished-isoform count
table. Bounded original-RNA replays test selected coding SNVs in both cell
lines; [HCC1395_ATLAS.md](HCC1395_ATLAS.md) and
[COLO829_ATLAS.md](COLO829_ATLAS.md) extend original-read status to every
selected small-variant source record and screen their SV sources.

## Deliverables

- [COLO829 whole-run atlas](COLO829_ATLAS.md): every SMaHT small variant has
  original-RNA coverage and allele status, ranked local protein contexts are
  listed separately, and every SV benchmark record has a breakpoint evidence
  screen. Its 10 MB indexed source-read BAM is committed under
  `tests/data/colo829_whole_run/` for offline replay.
- [HCC1395 whole-run atlas](HCC1395_ATLAS.md): every SEQC2 SNV/indel has an
  original direct-RNA allele status and local protein context outcome; every
  Talsania SV has a breakend screen. A bounded source-read BAM under
  `tests/data/hcc1395_whole_run/` replays informative calls offline.
- `somatic_call_inventory.tsv.gz`: 89,112 source records, preserving source
  identity, class, coordinates, allele where available, DNA filter/support,
  and published HCC1395 fusion annotation where it directly matches the SV.
  Its `rna_status=not_evaluated` column refers to the original inventory date;
  use `original_rna_protein_calls.tsv` for the bounded subset and the
  whole-run atlas tables for each cell line.
- `colo829_polished_isoform_counts.tsv.gz`: all 68,553 transcript IDs from
  TAGET's `final.mapped_fl_count.txt`, their COLO829T/COLO829BL full-length
  counts, and the tumor GFF interval where present. The 51,387 tumor-positive
  IDs all match a tumor GFF transcript. `sequence_allele_status` is
  `not_evaluated`; these counts are isoform-level full-length counts, not
  allele-spanning read counts or Isovar fragment support.
- `colo829_polished_snv_alleles.tsv.gz`: one row per callable SMaHT SNV and
  TAGET representative sequence. `allele=1` means the polished base matches
  the alternate allele; `0` means it matches the reference allele; `?` means
  it matches neither. Only primary alignments with mapping quality at least 20
  and a base covering the position appear. Tumor and normal full-length cluster
  counts remain separate columns.
- `colo829_polished_snv_summary.tsv.gz`: all 44,005 selected SMaHT SNVs,
  including those with no callable polished isoform. `not_callable` means no
  qualifying representative covered the site; it never means reference.
- `build_inventories.py` beside this report rebuilds the source-call and
  polished-count inventories from the downloaded files. It retains source
  records rather than merging independently defined benchmarks.
- `build_polished_alleles.py` rebuilds the two SNV allele tables from the
  versioned reference, SMaHT VCF, TAGET count table, and spliced BAM.
- `original_rna_protein_calls.tsv`: selected SEQC2 and SMaHT coding SNVs,
  separate original RNA read counts and Isovar protein calls per run, plus
  separate TAGET polished-sequence allele counts where available. `ref_reads`,
  `alt_reads`, and `other_reads` count Isovar allele-read objects after its
  default mate handling; `alt_source_reads` sums their original read counts.
  `protein_supporting_fragments` is support for the reported translation.
- `original_rna_alt_witnesses.tsv`: original alternate-read IDs and their
  source-read counts, linked to the indexed BAMs under
  `tests/data/hcc1395_colo829_original/`.
- `build_original_rna_report.py` rebuilds the selected-call and witness tables
  offline from the bounded BAMs and Ensembl 87 annotation.
- `extract_original_rna_vcfs.py` rebuilds the two bounded VCF fixtures from
  the exact SEQC2 and SMaHT source records. SEQC2 uses the compound filter
  `PASS;HighConf`; Varcode requires `only_passing=False` to retain those
  source records when loading the test fixture.
- `colo829_original_spot_ids.tsv` pins the 2,067 original SRA pair IDs used
  for the COLO829 fixture; `extract_colo829_original_pairs.py` recovers those
  pairs using SRA Toolkit `fastq-dump`.

To rebuild, install `openpyxl`, download the linked files below, and place
them in one directory under these names:

| Local name | Source file |
|---|---|
| `hcc1395-seqc2-highconf-snv.vcf.gz` | SEQC2 SNV VCF |
| `hcc1395-seqc2-highconf-indel.vcf.gz` | SEQC2 indel VCF |
| `13059_2022_2816_MOESM4_ESM.xlsx` | Talsania Table S4 |
| `13059_2022_2816_MOESM9_ESM.xlsx` | Talsania Table S9 |
| `Truthset_SNV_COLO829.vcf.gz` | SMaHT SNV VCF |
| `Truthset_Indel_COLO829.vcf.gz` | SMaHT indel VCF |
| `colo829-sv-grch38.vcf` | Paulin Paz 2024 GRCh38 v2 VCF |
| `colo829-sv-2022-hg38lifted.vcf` | Valle-Inclán 2022 lifted VCF |
| `colo829-fl-count.txt` | TAGET `final.mapped_fl_count.txt` |
| `colo829-tumor.gff` | TAGET `final.COLO829T.gff` |
| `colo829-final.rep.fq` | TAGET `final.rep.fq` |
| `GRCh38-no-alt.fa.gz` | NCBI GRCh38 no-alt analysis set (linked below) |

Run `python reports/hcc1395_colo829/build_inventories.py --input-dir
/path/to/downloads` from the repository root. Compare the source SHA-256
values below before interpreting the output.

## Source coverage

| Sample | Source | Records in inventory | Scope |
|---|---|---:|---|
| HCC1395 | [SEQC2 v1.2.1 SNV](https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/seqc/Somatic_Mutation_WG/release/latest/high-confidence_sSNV_in_HC_regions_v1.2.1.vcf.gz) | 39,447 | 37,330 `PASS;HighConf`; 2,117 `PASS;MedConf`; GRCh38.d1.vd1, high-confidence regions |
| HCC1395 | [SEQC2 v1.2.1 indel](https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/seqc/Somatic_Mutation_WG/release/latest/high-confidence_sINDEL_in_HC_regions_v1.2.1.vcf.gz) | 1,625 | 1,467 `PASS;HighConf`; 158 `PASS;MedConf` |
| HCC1395 | [Talsania et al. 2022 Table S4](https://link.springer.com/article/10.1186/s13059-022-02816-6) | 1,788 | Consensus SVs: 717 DEL, 551 INS, 230 DUP, 146 TRA, 133 INV, 11 BND. Table S9 names 24 RNA-seq fusion associated SVs; 20 IDs match Table S4, 4 do not. The table does not provide per-isoform sequences or read counts. |
| COLO829 | [SMaHT COLO829/BL SNV](https://github.com/parklab/SMaHT_SNV_COLO829BLT50_HAPMAP/tree/main/Resource/COLO829BLT50_reference_set) | 44,005 | Tumor-specific, matched-normal reference; `PASS`; GRCh38 no-alt |
| COLO829 | [SMaHT COLO829/BL indel](https://github.com/parklab/SMaHT_SNV_COLO829BLT50_HAPMAP/tree/main/Resource/COLO829BLT50_reference_set) | 2,059 | Same source and build |
| COLO829 | [Paulin Paz et al. SV v2](https://zenodo.org/records/13917379) | 54 | GRCh38 SV benchmark records, selected from replicates |
| COLO829 | [Valle-Inclán et al. 2022 SV](https://zenodo.org/records/4716169) | 134 VCF records | Earlier 68-event truth set lifted to GRCh38; breakends and other events can occupy multiple VCF records. Keep separate from v2. |
| COLO829 | [TAGET polished isoforms](https://zenodo.org/records/8319497) | 68,553 transcript-count rows | 51,387 have COLO829T full-length count > 0; 56,568 have COLO829BL count > 0. RNA source is a different study/culture from the SMaHT DNA truth set. |

## Original RNA replay and scope

The [HCC1395 ONT direct RNA run SRR23456709](https://www.ncbi.nlm.nih.gov/sra/SRR23456709)
has 1,049,452 reads. We extracted its first 100,000 original reads with NCBI
SRA Toolkit 3.4.1, aligned them to the NCBI GRCh38 no-alt analysis set with
minimap2 2.26 (`-ax splice -uf -k14`), and retained alignments overlapping the
OSTC and PRDX5 loci in an indexed BAM. Both SEQC2 variants have alternate
reads: OSTC F9L has 9, PRDX5 F157L has 18 under Isovar's default read rules.
The matched-normal SEQC2 VCF records report tumor/normal VAF of 0.490/0.001
for OSTC and 0.442/0.002 for PRDX5; these are DNA call-set values, not RNA
fractions.
Default Isovar translation produces no protein sequence from these direct RNA
subsets. A one-read, relaxed-coverage translation is possible but is not
reported as a supported mutant protein call. This is a measured limitation of
these noisy long-read subsets, not absence of the alleles.

We also streamed the first 100,000 original reads of the independent
[HCC1395 ONT direct cDNA run SRR23456710](https://www.ncbi.nlm.nih.gov/sra/SRR23456710)
through minimap2 (`-ax splice -k14`) and retained the same loci. At SEQC2
PRDX5 F157L, Isovar finds 41 reference, 24 alternate, and 2 other reads;
the top mutant protein is `GKVRLLADPTGALGKETDLLLDDSL`, with leucine at
zero-based `[12,13)` and six supporting fragments. OSTC F9L has two alternate
cDNA reads but no default translation. The Talsania Table S9 fusion names
remain external RNA evidence for 24 events, not original-read Isovar
mutant-isoform support or a complete protein ORF.
An [independent HCC1395 proteogenomics study](https://pmc.ncbi.nlm.nih.gov/articles/PMC11419083/)
reports PRDX5 F-to-L and OSTC F-to-L variant peptides from SEQC2 truth calls.
That peptide observation supports the mutated residue but does not validate
the entire Isovar protein context or establish that these RNA cultures are
identical to the study's proteomics culture.

The [COLO829 original paired-end RNA run SRR8615617](https://www.ncbi.nlm.nih.gov/sra/SRR8615617)
is a separate CCLE culture from both the SMaHT tumor/normal DNA and TAGET
long-read RNA. Its archived alignment uses **GRCh37**, while the SMaHT calls
use **GRCh38**. We mapped eight selected SMaHT coding sites through the
[UCSC hg38ToHg19 chain](https://hgdownload.soe.ucsc.edu/goldenPath/hg38/liftOver/hg38ToHg19.over.chain.gz),
queried the SRA `PRIMARY_ALIGNMENT` table at those positions, resolved each
alignment's `SEQ_SPOT_ID`, recovered 2,067 original read pairs, and realigned
them with minimap2 (`-ax sr`) to the NCBI GRCh38 reference. Isovar translates
seven of the eight selected alleles from those original reads. The extraction
depends on the archive's GRCh37 alignments, so its counts are bounded evidence
from those mapped spots, not exhaustive whole-run coverage.
The SMaHT VCF marks these records `MUT_ORIGIN=COLO829T` relative to its
COLO829BL matched normal; the TAGET normal isoform counts come from a distinct
library and cannot redefine that DNA baseline.
The spot table regenerates paired FASTQs with SHA-256 values
`1c5d4dc8795f1ff129367485a7dc39c405f4f8341b28f31bdf3a49557351b308`
(R1) and `914b2cd16a128b7255720f12349133a3efe5eef61988950d178552cfa7bc9fa0`
(R2).

For example, COLO829 NAXD `chr13:110634738:C>T` (Ensembl 87 p.P205S) has
104 reference and 81 alternate Isovar reads, representing 94 original source
reads. Isovar's top sequence is `LWLVAQQPALIHGYRKAVLTSNHVEFSRLYDAV`, with
serine at the zero-based interval `[20,21)` and 59 supporting fragments. The
independent TAGET comparison has four tumor-positive alternate representative
isoforms and zero normal alternate representatives at this site. These are
different assays and counts. SMG5, TNPO1, and ARHGAP5 also have alternate
TAGET normal representatives despite being selected from a SMaHT somatic DNA
set; the distinct cell-line cultures and source definitions must remain visible.

The bounded BAMs and `original_rna_candidates.tsv` are the offline test data.
Run `python reports/hcc1395_colo829/build_original_rna_report.py` from the
repository root to regenerate the two original RNA tables. The protein calls
use Ensembl 87, Isovar's default read and translation settings, and zero-based
mutation intervals. An `NA` protein sequence with `no_translation` is an
evaluated, allele-supported outcome. The full DNA inventory and polished
isoform files retain their original `not_evaluated` values; the newer atlases
provide small-variant outcomes and separate SV screens. Neither atlas asserts
an SV-derived protein.
Protein positions depend on the selected reference transcript: for example,
Varcode's top SYT14 effect is p.S682F under Ensembl 87 and p.S354F under
Ensembl 95 for the same genomic call.
To regenerate the VCFs from the pinned source downloads, run
`python reports/hcc1395_colo829/extract_original_rna_vcfs.py --hcc-vcf
hcc1395-seqc2-highconf-snv.vcf.gz --colo-vcf
Truthset_SNV_COLO829.vcf.gz`.

With NCBI SRA Toolkit 3.4.1, minimap2 2.26, samtools 1.21, and the NCBI
GRCh38 no-alt FASTA, regenerate the original RNA alignments as follows
(`$REF` is the FASTA path, `$DUMP` is the `fastq-dump` executable, and
`$OUT` is a scratch directory):

```sh
python reports/hcc1395_colo829/extract_colo829_original_pairs.py \
  --fastq-dump "$DUMP" --output-dir "$OUT/colo829"
minimap2 -t 4 -ax sr "$REF" "$OUT/colo829/R1.fastq" "$OUT/colo829/R2.fastq" \
  | samtools sort -o "$OUT/colo829/original.grch38.bam"
samtools index "$OUT/colo829/original.grch38.bam"

"$DUMP" -N 1 -X 100000 -Z SRR23456709 \
  | minimap2 -t 4 -ax splice -uf -k14 "$REF" - \
  | samtools view -h -b -L reports/hcc1395_colo829/hcc1395_candidate_regions.bed \
      -o "$OUT/hcc1395_direct_unsorted.bam"
samtools sort -o "$OUT/hcc1395_srr23456709.bam" "$OUT/hcc1395_direct_unsorted.bam"
samtools index "$OUT/hcc1395_srr23456709.bam"

"$DUMP" -N 1 -X 100000 -Z SRR23456710 \
  | minimap2 -t 4 -ax splice -k14 "$REF" - \
  | samtools view -h -b -L reports/hcc1395_colo829/hcc1395_candidate_regions.bed \
      -o "$OUT/hcc1395_cdna_unsorted.bam"
samtools sort -o "$OUT/hcc1395_srr23456710.bam" "$OUT/hcc1395_cdna_unsorted.bam"
samtools index "$OUT/hcc1395_srr23456710.bam"
```

TAGET supplies full-length counts, a tumor GFF, and polished representative
transcript sequences. A preliminary length audit found only 220 of the first
702 GFF-matched representative sequences had length equal to the sum of their
GFF exons. Thus a direct genomic-coordinate to polished-sequence offset is
invalid for most isoforms. The spliced alignment below tests **SNV alleles in
polished sequences**. Even when an allele is present, the full-length count is
support for the isoform cluster, not a count of original reads crossing that
allele.

## Binary allele observations from polished COLO829 isoforms

We aligned all 68,553 TAGET `final.rep.fq` representatives to the [NCBI
GRCh38 no-alt analysis set](https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.15_GRCh38/seqs_for_alignment_pipelines.ucsc_ids/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gz)
with minimap2 2.26 (`-ax splice:hq -uf`) and samtools 1.21, then checked all
44,005 SMaHT COLO829 SNVs against the aligned sequence. The VCF reference base
was verified against the FASTA at every site. A base is callable only on a
primary alignment with mapping quality at least 20. An uncovered, skipped,
deleted, or ambiguously mapped site has no 0/1 row.

| Observation | Count |
|---|---:|
| SNVs with any callable polished isoform | 425 |
| SNVs with a callable tumor-positive isoform | 330 |
| SNVs with a tumor-positive alternate isoform | 132 |
| Callable SNV/isoform pairs | 1,642 |
| Pairs with alternate (1), reference (0), other (?) base | 263 / 1,378 / 1 |
| Tumor-positive pairs with alternate (1), reference (0), other (?) base | 263 / 880 / 1 |

For example, at SMaHT `chr1:948136:G>A`, TAGET representative `PB.7.4` has
reference `G` (`0`), while `PB.7.5` and `PB.7.8` have alternate `A` (`1`). Their
tumor full-length cluster counts are 68, 2, and 7, respectively. A 5.1 KB
source-derived BAM with these three alignments is the regression fixture in
`tests/data/colo829_polished/`.

These are **sequence observations**, not allele-spanning original-read counts,
mutant protein calls, or a statement that the same passages carry a mutation.
The SMaHT DNA and TAGET RNA come from different studies. In particular, 125
alternate-sequence pairs also have a nonzero TAGET normal cluster count; this
audit cannot explain those shared counts. Indels and SV breakpoints still need
separate alignment and haplotype checks. The selected HCC1395 and COLO829
original RNA results are in the source-separated table above.

To rebuild the sequence audit after downloading the files listed here and the
NCBI FASTA, run:

```sh
minimap2 -ax splice:hq -uf -t 6 GRCh38-no-alt.fa colo829-final.rep.fq \
  | samtools sort -o colo829-polished.grch38.bam -
samtools faidx GRCh38-no-alt.fa
python reports/hcc1395_colo829/build_polished_alleles.py \
  --bam colo829-polished.grch38.bam \
  --vcf Truthset_SNV_COLO829.vcf.gz \
  --reference GRCh38-no-alt.fa \
  --counts colo829-fl-count.txt
```

The selected DNA call sets are bounded by their own reference regions, filters,
variant representations, and cell passages. They are not a claim that all
biological mutations in either cell line are known. A negative RNA result will
require callable coverage; zero supporting reads alone is not evidence of
absence. SVs need breakpoint-spanning RNA analysis and their transcript/ORF
interpretations cannot be inferred from DNA coordinates alone.

## Remaining analysis

The HCC1395 and COLO829 whole-run atlases assign RNA support statuses to
their small variants and screen their source SV records for breakend evidence.
SVs still require path-aware breakpoint RNA analysis and separate ORF
reconstruction. A full-length mutant isoform also
requires read linkage across its splice chain; the local protein fragments in
this report do not assert that linkage. Direct RNA and cDNA disagreement at
translation stage should be evaluated with read-quality-aware consensus, not
silently filled from DNA annotation or polished representative counts.

## Input SHA-256

| Input | SHA-256 |
|---|---|
| HCC1395 SEQC2 SNV gzip | `0c49f7e1519fc7860b773f79bdcf07453fe4def8303707f34542471d92494c9d` |
| HCC1395 SEQC2 indel gzip | `db7dc8be2ec04bab830ff1411146c7d538c217ff684e874b85d98a3dbbda91cf` |
| HCC1395 SV Table S4 xlsx | `11e8d6393b83f0e0b7fdb0b13a57376663fc9662116d78c89a0d5fa530df096b` |
| HCC1395 RNA fusion Table S9 xlsx | `f061444e4fefdc36894553ba99f46bfcff6fa04d5957e48b69f6fd3647968043` |
| COLO829 SMaHT SNV gzip | `e57f77c8bc280bbad2255856c4062348f92f25da38726f367322c52736578f31` |
| COLO829 SMaHT indel gzip | `e8d0fb67925a2365b485cd65feab8dfdc12bbe54c0bc0bdec44c1846b518b685` |
| COLO829 SV v2 GRCh38 VCF | `8ad075bf43aa71d89c591396569e813959ada44a30f171982c2d299411290bac` |
| COLO829 earlier SV GRCh38 lifted VCF | `a8f38454367cf8837a7dba1d29c9bc338a8e62cadfee6d2489931f8105f0dcdc` |
| COLO829 TAGET full-length count table | `c54e284479ddb78714a30e03dfdec04e38b2d45fbb26f0a96717d83e5e9d80d9` |
| COLO829 TAGET tumor GFF | `506488c868150e35358514fd1b9528ca4070ed0f4297e50d42364ad6594f1b94` |
| COLO829 TAGET representative FASTQ | `6b008d7610c03bdcc62e716bd365f4b448badc9c77554c58574068c7c9ab487b` |
| NCBI GRCh38 no-alt FASTA gzip | `fb4243ebb014caf27111f24dd62b7ce42160f28581da6f8fcd6cba5977778d02` |

## Original RNA provenance

The [ENA run file report](https://www.ebi.ac.uk/ena/portal/api/filereport?accession=SRR8615617&result=read_run&fields=run_accession,fastq_ftp,fastq_bytes,fastq_md5,read_count,base_count&format=tsv)
lists SRR8615617 paired FASTQ MD5 values
`799b4d8f1a2743df0d74aa34936253d3` and
`41c9d6bcf4c330655ccb0cb33a7e260d`. The
[HCC1395 ENA report](https://www.ebi.ac.uk/ena/portal/api/filereport?accession=SRR23456709&result=read_run&fields=run_accession,fastq_ftp,fastq_bytes,fastq_md5,read_count,base_count&format=tsv)
lists SRR23456709 FASTQ MD5 `1f5d2a60c3158ae12e573a316bd846bb`.
The [cDNA ENA report](https://www.ebi.ac.uk/ena/portal/api/filereport?accession=SRR23456710&result=read_run&fields=run_accession,fastq_ftp,fastq_bytes,fastq_md5,read_count,base_count&format=tsv)
lists SRR23456710 MD5 `223d07e4e91a63a6380349758aa3ff09`.
These are archive metadata for the full FASTQs; we did not download or
rehash those full files. The bounded BAMs committed here have SHA-256 values:

| Indexed original-read BAM | SHA-256 |
|---|---|
| HCC1395 direct RNA SRR23456709 | `b5d6b883b0a1982d3f53d582b14ff95656232f9f9ece626510c4cd899b3b9a09` |
| HCC1395 direct cDNA SRR23456710 | `d46cd11d4da900634814c423345374210495f4417cf65466973c8a04b0081bb4` |
| COLO829 paired RNA SRR8615617 | `459afe345fdae603123577774836399f64cc40ba5cfcdae78612c7f7f7643e55` |

NCBI's [public SRA data policy](https://www.ncbi.nlm.nih.gov/sra/docs/sequence-data-processing/)
states that public sequence data may be retrieved and redistributed. The
bounded BAMs retain only source-derived reads around the selected loci.
