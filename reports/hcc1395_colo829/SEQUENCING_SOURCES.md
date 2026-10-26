# Sequencing, nomination and reference sources

This index covers the Sid, COLO829 and HCC1395 report set linked from
[THREE_SAMPLE_SUMMARY.md](THREE_SAMPLE_SUMMARY.md), including the focused
[CDKN2A](focused_proteins/COLO829_CDKN2A_frameshift.md) and
[EIF3K-CYP39A1](focused_proteins/HCC1395_EIF3K_CYP39A1_fusion.md) alignments.
Source libraries, timepoints and cultures remain separate. Derived BAMs,
selected windows, polished representatives and multiple reports from one
library are not additional biological replicates.

## Original RNA libraries

| Sample / library | Original source identifier | Used in these reports | Scope and provenance |
|---|---|---|---|
| COLO829 Illumina paired cDNA RNA | [SRR8615617](https://www.ncbi.nlm.nih.gov/sra/SRR8615617), PRJNA523380, SAMN10988208 | Whole-run atlas, CDKN2A/FHL2 windows, SV screen and unified SV audit | 115,269,780 pairs / 230,539,560 101-base reads. Complete source run aligned before selected-site extraction. CCLE RNA culture differs from the benchmark DNA and TAGET cultures. [Run metadata](https://www.ebi.ac.uk/ena/portal/api/filereport?accession=SRR8615617&result=read_run&fields=run_accession,sample_accession,study_accession,library_strategy,library_source,library_layout,read_count,base_count,fastq_bytes,fastq_md5&format=tsv). |
| HCC1395 ONT direct RNA | [SRR23456709](https://www.ncbi.nlm.nih.gov/sra/SRR23456709), PRJNA934933, SAMN33292477 | Whole-run atlas, selected SV replay, independent EIF3K-CYP39A1 evidence | 1,049,452 original reads. Whole-run atlas follows complete-run alignment. The fusion check separately retains 21 original candidates and 42 local alignments from the same library; these are not additional replicates. [Run metadata](https://www.ebi.ac.uk/ena/portal/api/filereport?accession=SRR23456709&result=read_run&fields=run_accession,sample_accession,study_accession,library_strategy,library_source,library_layout,read_count,base_count,fastq_bytes,fastq_md5&format=tsv). |
| HCC1395 ONT direct cDNA | [SRR23456710](https://www.ncbi.nlm.nih.gov/sra/SRR23456710), PRJNA934933 | Earlier bounded OSTC/PRDX5 original-RNA comparison | First 100,000 reads streamed; selected loci retained. This library is separate from direct RNA and does not enter its whole-run funnel. [Methods and results](README.md#original-rna-replay-and-scope). |
| HCC1395 Illumina tumor RNA | [INTEGRATE public example](https://sourceforge.net/p/integrate-fusion/wiki/example/), `test-data.tar.gz`: `unmapped.bam` and `accepted_hits.bam` | EIF3K-CYP39A1 coding reconstruction | Archive-specific source ID; no SRA accession is asserted for this subset. 185 original fusion-enriched unmapped records screened, plus 19 original mapped coding-tiling records. Junction support is not whole-library frequency. [Fixture provenance](../../tests/data/hcc1395_eif3k_cyp39a1/README.md), [primary paper](https://pmc.ncbi.nlm.nih.gov/articles/PMC4691743/). |
| Sid bulk tumor RNA T2 | Osteosarc product `40a3ceff825d…`; `rna-seq/reprocessed/SARC0277/SARC0277.Aligned.sortedByCoord.out.bam` | Selected indel/frameshift and SV/ORF reports | SARC0277 T2. Regional acquisition and recorded caps; source keys and full product IDs in the [input index](../../notes/data/sid-neoorf/input-files.tsv). |
| Sid tagged ONT tumor RNA T2 | Osteosarc product `58f81602ec2f…`; `IPISRC044_T2_sclrs_ONT.tagged.bam` | Selected indel/frameshift and SV/ORF reports | Original registered key is `ONT/IPISRC044_ONT_upload/IPISRC044_ONT/processed/IPISRC044_T2_sclrs_ONT/IPISRC044_T2_sclrs_ONT/IPISRC044_T2_sclrs_ONT.tagged.bam`. Regional inputs can have incomplete mate/supplementary recovery. |
| Sid Kamil 10x tumor RNA T2 | Osteosarc product `6e9147690525…`; `kamil/tumor/output/T2/outs/per_sample_outs/T2/sample_alignments.bam` | Selected indel/frameshift and SV/ORF reports | Regional acquisition; read fragments and cellular labels remain different measurements. Repeated processing products do not constitute independent libraries. |
| Sid PacBio tumor RNA T1 | Osteosarc product `8275c5f21164…`; `pacbio/IPISRC044_T1_sclrs_live_pbmm2_mapped.bam` | Selected indel/frameshift and SV/ORF reports | Distinct timepoint from T2. Regional recovery can be incomplete; some original query records lack qualities, which remain explicitly unassessed. |

Sid source acquisition is owned by [Osteosarc](https://github.com/iskandr/osteosarc),
using its registered source inventory and recorded acquisition receipts. No SRA
accession is invented for those product keys. The [selected-read methods](../../notes/sid-neoorf-read-subsets.md),
[four-platform pilot](../../notes/sid-sv-four-platform-pilot.md),
[source-date checks](../../notes/data/sid-sv-recovery/source-date-checks.json),
and [bundle file hashes](../../notes/data/sid-neoorf/bundle-files.tsv) record
identity, selection and completeness. The full catalogue of 77 candidate
processing products has not been audited; the report uses the four named
products. The [report download manifest](../../notes/data/sid-neoorf-event-report/downloads.json)
pins its sequence/evidence ledger; it is a derived report, not another RNA library.

## DNA nominations and matched-normal evidence

| Sample | Source | Contribution and limits |
|---|---|---|
| HCC1395 / HCC1395BL | SEQC2 v1.2.1 [SNV VCF](https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/seqc/Somatic_Mutation_WG/release/latest/high-confidence_sSNV_in_HC_regions_v1.2.1.vcf.gz) and [indel VCF](https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/seqc/Somatic_Mutation_WG/release/latest/high-confidence_sINDEL_in_HC_regions_v1.2.1.vcf.gz) | 39,447 SNVs and 1,625 indels, GRCh38.d1.vd1. HighConf and MedConf source labels retained. VCF tumor/normal VAF is DNA evidence, not an RNA fraction. [Primary study](https://pmc.ncbi.nlm.nih.gov/articles/PMC8532138/). |
| HCC1395 / HCC1395BL | [Talsania et al. 2022](https://link.springer.com/article/10.1186/s13059-022-02816-6), Tables S4 and S9 | 1,788 consensus SV records and published RNA-fusion annotations. These annotations are not original-read reconstruction evidence. |
| HCC1395 / matched normal | [INTEGRATE example archive](https://sourceforge.net/p/integrate-fusion/wiki/example/) tumor and normal genomic BAMs | Original alignments from GRCh37 chr19 `[39124000,39125000)` and chr6 `[46608000,46609000)` retained in separate fixtures. EIF3K-CYP39A1: 18 tumor cross-gene fragments / 0 normal. Not a base-resolved DNA breakpoint. [Exact fixture hashes](../../tests/data/hcc1395_eif3k_cyp39a1/README.md). |
| COLO829 / COLO829BL | [SMaHT COLO829BLT50 reference set](https://github.com/parklab/SMaHT_SNV_COLO829BLT50_HAPMAP/tree/main/Resource/COLO829BLT50_reference_set), `Truthset_SNV_COLO829.vcf.gz`, `Truthset_Indel_COLO829.vcf.gz` | 44,005 SNVs and 2,059 indels; PASS, COLO829T-origin nominations relative to COLO829BL; GRCh38 no-alt. CDKN2A comes from the indel set. |
| COLO829 | [Paulin Paz 2024 GRCh38 v2](https://zenodo.org/records/13917379) | 54 SV VCF records, preserved separately from the older benchmark. |
| COLO829 | [Valle-Inclan 2022](https://zenodo.org/records/4716169), lifted VCF | 134 VCF records from the earlier 68-event set. Paired BNDs are separate source rows; these are not 134 independent events. |
| Sid | Frozen vaccine/small-variant metadata and original UCLA/Personalis somatic calls, including T0, T1, T2 and T1 organoid | The [55 event sheets](../../notes/data/sid-neoorf-event-report/README.md) retain per-event original records, timepoint-specific DNA depth/VAF or VF, and source discrepancies. RNA-only caller nominations remain labelled RNA-only; missing DNA origin is not inferred from an RNA join. The supplemental ATP5MG-KMT2A nomination is pinned to its [original-RNA input fixture](../../tests/data/fusions/coding-corpus/ATP5MG--KMT2A.input.json.gz). |

## Polished and unavailable RNA sources

| Source | Use / availability |
|---|---|
| [TAGET COLO829 dataset](https://zenodo.org/records/8319497): `final.rep.fq`, `final.COLO829T.gff`, `final.mapped_fl_count.txt` | Polished representatives and separate COLO829T/COLO829BL full-length cluster counts. 68,553 count rows. Used in the independent polished-allele and SV3 cross-checks; never counted as original RNA molecules or pooled with SRR8615617. [Cross-check](colo829_sv3_polished_check.md). |
| Original 2019 COLO829 tumor and normal Iso-Seq | Declared as two separate unavailable sources in the [unified manifest](unified_sv_audit/manifest.json). The original host timed out; alternate host returned 404 on 2026-10-26. [Recovery #574](https://github.com/openvax/isovar/issues/574). These rows are unassessed; the polished archive does not supply the missing originals. |

## Reference sequences and derived inputs

- **DDX17 follow-up:** [pinned original-RNA audit](ddx17_audit/README.md) uses
  SRR8615617 only: all 603 original site-selected pairs, plus regional source-index
  recovery, yield 27,309 distinct pairs. Original SRA GRCh37 placements and fresh
  STAR 2.7.11b GRCh38 local placements are kept distinct. The 50-kb GRCh38 genomic
  slice, Ensembl 87 transcript sequences, original FASTQs, alternative placements,
  selection bounds and SHA-256 receipts are in `ddx17_audit/manifest.json`.
  This bounded investigation is not a new culture or a whole-library splice audit.

- **Cell-line genome:** NCBI GRCh38 no-alt analysis FASTA; uncompressed SHA-256
  `9cce8b926416dd96b152deea85188495b75f7ac8d634cc723a017067be8702b7`.
  Download URL and source archive hashes are in the [inventory README](README.md#input-sha-256).
- **Cell-line transcripts and proteins:** Ensembl release 87, GRCh38. The
  [reference-window table](example_reference_contexts.tsv) pins each comparison
  transcript, residue offset, sequence and reference-protein SHA-256. The
  protein FASTA used for the reference-self screen has SHA-256
  `d873c199d5dbf20869db3f2cbd6ac726259215d91c2782f9f1ddcad4bb117cf9`.
- **HCC fusion:** GRCh37-lite exon/cDNA models supplied by the INTEGRATE archive.
  Ensembl 87 CDS offsets transfer only for identical whole-transcript sequences.
  ONT alignments use separate GRCh38 exon models; coordinates are never mixed.
  All eight inputs have hashes in the focused report JSON and fixture README.
- **Cell-line ncORF annotation:** orfdata 0.34.1, human GRCh38 catalogue
  2026-10-11.1, 7,264 reference ncORF models based on Ensembl 101. Catalogue
  SHA-256 `06bc0bf3f81c25a7356f71a38581ff5704c05dcb75752454bc1ca75075ab8c09`.
  These are coding-block overlaps, not reconstructed mutant ncORF proteins.
- **Sid references:** pinned Ensembl 115 GRCh38 transcripts/proteins. The
  [Sid report methods](../../notes/sid-neoorf-event-report.md) and its input
  manifests retain reference identities and exact/I-L-collapsed comparison
  rules. Sid and cell-line novelty counts use different reference releases.
- **Coordinate conversion:** the early COLO829 spot-selected replay uses the
  [UCSC hg38ToHg19 chain](https://hgdownload.soe.ucsc.edu/goldenPath/hg38/liftOver/hg38ToHg19.over.chain.gz)
  to find original SRA spots before GRCh38 realignment. The later complete-run
  atlas does not use that older alignment to preselect reads.
- **Reproducible derivatives:** [whole-run COLO829 fixture](../../tests/data/colo829_whole_run/README.md),
  [whole-run HCC1395 atlas](HCC1395_ATLAS.md),
  [original selected-call fixtures](../../tests/data/hcc1395_colo829_original/README.md),
  [fusion fixture](../../tests/data/hcc1395_eif3k_cyp39a1/README.md), and
  [unified SV manifest](unified_sv_audit/manifest.json) identify every retained
  source, reference file, selection and budget. The SV snapshots have individual
  SHA-256 receipts. None of these derivatives is an independent library.

RNA sequence support does not measure translation. Missing inputs, incomplete
acquisition and bounded discovery remain unknown rather than negative evidence.
