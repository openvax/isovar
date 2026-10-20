# HCC1395 and COLO829: source-call and RNA-isoform audit

Date: 2026-10-20. This is a **source inventory and COLO829 polished-sequence
SNV audit**, separate from the Isovar WLOS_0001 test. It enumerates every
record in the selected public DNA call sets and every row in TAGET's COLO829
polished-isoform count table. No original HCC1395 or COLO829 RNA read set has
been replayed against these calls, and no variant-specific mutant protein
sequence is asserted.

## Deliverables

- `somatic_call_inventory.tsv.gz`: 89,112 source records, preserving source
  identity, class, coordinates, allele where available, DNA filter/support,
  and published HCC1395 fusion annotation where it directly matches the SV.
  All `rna_status` values are `not_evaluated`.
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

## What the RNA evidence presently says

- **HCC1395:** Original long-read tumor RNA is available as
  [PRJNA934933](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA934933), including
  ONT direct RNA SRR23456709. Its ENA FASTQ is 1,167,327,234 compressed bytes.
  An SRA `sam-dump --aligned-region` trial emitted unaligned whole-run data,
  so there is no verified region-only aligned fixture from that accession.
  The Talsania Table S9 fusion names are external RNA evidence for 24 events,
  not original-read Isovar mutant-isoform support or a complete protein ORF.
- **COLO829:** TAGET supplies full-length counts, a tumor GFF, and 450 MB of
  polished representative transcript sequences. A preliminary length audit
  found only 220 of the first 702 GFF-matched representative sequences had
  length equal to the sum of their GFF exons. Thus a direct genomic-coordinate
  to polished-sequence offset is invalid for most isoforms. The spliced
  alignment below now tests **SNV alleles in polished sequences**. Even when an
  allele is present, the full-length count is support for the isoform cluster,
  not a count of original reads crossing that allele.

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
separate alignment and haplotype checks. HCC1395 RNA remains unevaluated.

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

## Next analysis needed for the requested mutant-isoform report

1. Pin the exact tumor RNA library and GRCh38 build for each line; align its
   original long reads and retain primary and alternative placements.
2. Replay each small variant and SV breakpoint with Isovar, recording allele
   reads/fragments, coverage, full versus local protein status, and all
   alternative isoforms. Mark uncallable and unsupported separately.
3. For COLO829, compare the polished-sequence observations above with original
   reads. Never substitute full-length cluster counts for variant-spanning
   witnesses. Extend the allele audit to indels and validated SV breakpoints.
4. Join the outputs to this complete source-record inventory and retain call
   set, reference build, sample passage, and source hash per row.

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
