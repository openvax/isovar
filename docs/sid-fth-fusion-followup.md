# FTH1–FTL junction predictions and longitudinal fusion calls

For the DNA-supported recurrent event and ORF priorities, see the
[curated somatic fusion shortlist](sid-recurrent-somatic-fusions.md).

The reconstructed 179-aa FTH1(1–17)–FTL(14–175) RNA ORF is absent from the available public vaccine target/peptide records. No explicit historical rejection reason was found; absence from the final list does not show that it was considered and rejected.

## New peptide prediction screen

Installed NetMHCpan 4.2c evaluated all 34 junction-spanning 8–11-aa windows against the five documented classical class-I alleles (170 predictions). None passed the default EL rank <2% threshold. The closest was DSEAAVNSL / HLA-A*01:01, EL rank 2.049% (BA affinity 14,626 nM). This is a model-specific screen, not proof that no CD8 epitope exists.

Installed NetMHCIIpan 4.3e evaluated fourteen junction-spanning 15-mers against DRB1*03:01, DRB1*08:01 and DRB3*01:01 (42 predictions). Five overlapping DRB3*01:01 windows were strong and one was weak at the default EL thresholds. They share the core YHQDSEAAV; these are one overlapping region, not six independent epitopes. The strongest 15-mer, VRQNYHQDSEAAVNS, had EL rank 0.07% and predicted affinity 147 nM. DRB1 had no default-threshold hits; DP/DQ were not screened.

The core and 22 of the 34 class-I windows are absent from both supplied canonical parent proteins. Whole-proteome novelty, protein expression, antigen presentation and T-cell recognition remain untested. This screen covers this exact reconstructed isoform, not every FTH1–FTL splice product. These new predictions cannot explain historic vaccine selection without its candidate and rejection logs.

[Patient HLA summary](https://osteosarc.com/dragen/hla/); [NetMHCpan documentation](https://services.healthtech.dtu.dk/services/NetMHCpan-4.2/); [NetMHCIIpan documentation](https://services.healthtech.dtu.dk/services/NetMHCIIpan-4.3/). Raw local outputs, input and software hashes, and the complete prediction grid are pinned in [the prediction ledger](data/sid-sv-recovery/fth-junction-predictions.json).

## FTH1–FTL recurrence

In the preserved filtered CTAT long-read table, the gene pair is called only in T2 (45, 18 and 12 reads at three distinct junctions). T1 and T3 have no filtered call; T0 is outside this long-read table. These are missing calls, not demonstrated zero expression. The two complete witnesses for our exact 179-aa ORF remain T2-only and are not the same exported read identifiers as the CTAT calls.

Across patients, a primary lung-cancer study reported FTH1–FTL in 3/49 patients with detected fusion calls (among 54 patients screened), and prior appearances in CCLE and TCGA. That establishes recurrence of the gene pair, without establishing recurrence of our exact junction/protein or a somatic DNA mechanism. CTAT FusionAnnotator retrieves catalogue tags by gene-pair name. [Primary study](https://pmc.ncbi.nlm.nih.gov/articles/PMC8002233/); [annotation implementation](https://github.com/FusionAnnotator/FusionAnnotator/blob/master/FusionAnnotator).

## Shared exact long-read junctions

The filtered T1/T2/T3 table contains 16 strand-aware genomic junctions shared across time points, belonging to 13 gene pairs. Below, comma-separated counts represent separate junctions; blank entries mean no filtered call. They must not be interpreted as zero RNA support or summed independent molecules. Matching genomic junctions do not establish identical ORFs.

| Gene pair | T1 supporting reads | T2 supporting reads | T3 supporting reads |
| --- | ---: | ---: | ---: |
| ENSG00000266869–RGS6 | 5 | 15 | — |
| ENSG00000279880–BPTF | 30 | 41 | 57 |
| ENSG00000289495–TRAC | 6 | 25 | 40 |
| FRG1HP–ENSG00000291075 | 7 | 14 | — |
| HAND2-AS1–LINC02269 | 5 | 10 | — |
| MAF–ENSG00000261722 | 17 | 31 | — |
| MROCKI–LINC02541 | 38 | 57 | — |
| OAZ1–SCO2 | — | 11 | 13 |
| TALAM1–ENSG00000280614 | — | 31, 20 | 19, 12 |
| TALAM1–ENSG00000280800 | — | 31, 20 | 19, 12 |
| TALAM1–ENSG00000281181 | — | 31, 20 | 19, 11 |
| TALAM1–RPL23A | — | 46 | 22 |
| TPST1–CRCP | 27 | 22 | — |

Only ENSG00000279880–BPTF and ENSG00000289495–TRAC have an identical filtered junction called at all three time points. The latter is immune-receptor-associated. TALAM1 partner ambiguity and other recurrent RNA/library artifacts require review; recurrence alone is not a tumor-specific protein finding.

The broader caller-independent DNA/RNA probe highlights FOXO3–STRADA/CCDC47, PARD3B–CDKN2B-AS1/CDKN2B, the TPST1/CRCP-region event, and GABBR1–SLC29A1 as events with RNA support at T1/T2/T3 and corroborating DNA at earlier samples. These labels and genomic-probe counts are a different evidence layer from filtered CTAT calls and complete ORF witnesses. [Published fusion evidence](https://osteosarc.com/fusions/).

[The longitudinal ledger](data/sid-sv-recovery/fusion-longitudinal-comparison.json) preserves each junction, source hash and per-sample count. [The CTAT source page](https://osteosarc.com/ctat_lr_fusion/) documents its filtered-call definition.

## What is pending for T3

T3 RNA fusion calling is complete (24 filtered calls). The CTAT page still says “UCLA — pending mapping,” while preserved BAM metadata identifies the T3 ONT library as UCSF / April 2025. The clinical T3 resection is dated April 17, 2025 on the [sample page](https://osteosarc.com/scrnaseq/); an exact library-to-specimen receipt was not newly verified. No T3 oncoanalyser tumor/normal counterpart is available in the CTAT comparison, so short-read/DNA concordance is unavailable. Neither gap prevents using the existing T3 RNA calls.

The separate unfinished Isovar recovery views are IGKC in the T2 ONT library. They do not imply that T3 fusion calling is waiting on this recovery.
