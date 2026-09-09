# Full-region vaccine RNA/protein audit

44 vaccine variants × 164 RNA alignment products = 7216 rows.

These are products, **not independent samples**. Technical reprocessings and unresolved pooled/blood libraries are not combined.

Counts below count source/variant rows, never pooled biological read totals. Missing acquisition and interrupted computation are not zero support.

## Default-path outcomes

| Outcome | Source/variant rows |
|---|---:|
| acquisition_failed | 442 |
| audit_timeout | 5 |
| different_top_protein | 16 |
| expected_top_protein | 749 |
| missing_index | 220 |
| no_alt_reads | 3740 |
| no_callable_allele | 568 |
| no_genomic_coordinates | 352 |
| no_overlapping_alignment | 664 |
| no_rna_candidate | 382 |
| no_translation | 78 |

## Per-variant coverage

These columns describe the default-returned protein, not hidden uncapped alternatives. Each count is the number of RNA products with that result; one product can represent pooled donors or a duplicate processing.

| Variant | Products with counts | With alt | With default protein | Expected default | With mutant 25-mer |
|---|---:|---:|---:|---:|---:|
| ABCF2-chr7-151218156 | 141 | 10 | 0 | 0 | 0 |
| ADGRF5-chr6-46856758 | 141 | 0 | 0 | 0 | 0 |
| ATRX-chrX-77599465 | 141 | 24 | 10 | 10 | 10 |
| BMP1-chr8-22192067 | 141 | 28 | 25 | 25 | 25 |
| BTD-chr3-15645183 | 141 | 23 | 19 | 19 | 17 |
| CD109-chr6-73818405 | 141 | 23 | 16 | 16 | 16 |
| CDC40-chr6-110228838 | 141 | 22 | 20 | 20 | 20 |
| CENPL-chr1-173811223 | 141 | 17 | 15 | 15 | 13 |
| CUL9-chr6-43186384 | 141 | 16 | 2 | 2 | 2 |
| DIAPH1-chr5-141526090 | 141 | 21 | 9 | 9 | 9 |
| DYNC1H1-chr14-101980529 | 141 | 52 | 44 | 44 | 44 |
| DYNC1H1-chr14-102030200 | 141 | 10 | 6 | 6 | 6 |
| EPG5-chr18-45925851 | 141 | 12 | 8 | 8 | 8 |
| EXOC4-chr7-133274996 | 141 | 110 | 76 | 76 | 63 |
| F8-chrX-154930861 | 141 | 10 | 4 | 4 | 4 |
| FBXO33-chr14-39431722 | 141 | 61 | 15 | 15 | 8 |
| GLIS3-chr9-3856149 | 141 | 6 | 4 | 4 | 4 |
| H1_2-chr6-26055824 | 141 | 19 | 19 | 19 | 15 |
| KDM3B-chr5-138391529 | 141 | 24 | 12 | 12 | 12 |
| KIF1C-chr17-5007293 | 141 | 14 | 12 | 12 | 12 |
| KRT18-chr12-52951767 | 141 | 8 | 4 | 4 | 4 |
| MAP2-chr2-209694768 | 141 | 2 | 0 | 0 | 0 |
| MT_ND5-chrM-12994 | 139 | 126 | 101 | 101 | 66 |
| MYO15B-chr17-75589953 | 141 | 2 | 0 | 0 | 0 |
| MYO9A-chr15-71994488 | 141 | 8 | 0 | 0 | 0 |
| NAV2-chr11-20080142 | 141 | 16 | 12 | 12 | 12 |
| NME1-chr17-51154443 | 141 | 39 | 39 | 37 | 31 |
| NR2F2-chr15-96332299 | 141 | 1 | 1 | 0 | 1 |
| NTF3-chr12-5494381 | 141 | 9 | 9 | 0 | 4 |
| PIP5K1A-chr1-151242178 | 141 | 21 | 18 | 18 | 11 |
| PRICKLE2-chr3-64198863 | 141 | 6 | 6 | 6 | 6 |
| PTH1R-chr3-46898732 | 141 | 10 | 10 | 10 | 10 |
| RNF213-chr17-80327830 | 141 | 27 | 18 | 18 | 18 |
| ROBO2-chr3-77607853 | 141 | 34 | 25 | 24 | 23 |
| SLC25A12-chr2-171813470 | 141 | 38 | 33 | 33 | 29 |
| SMC5-chr9-70298024 | 141 | 91 | 24 | 24 | 24 |
| TECPR1-chr7-98241127 | 141 | 33 | 14 | 14 | 5 |
| TMRO-chr9-97910412 | 141 | 55 | 28 | 28 | 28 |
| VPS13B-chr8-99861923 | 141 | 25 | 5 | 5 | 5 |
| VPS72-chr1-151185547 | 141 | 45 | 27 | 27 | 25 |
| WWC1-chr5-168399513 | 141 | 17 | 7 | 7 | 7 |
| ZNF436-chr1-23362762 | 141 | 30 | 27 | 27 | 19 |
| ZNF674-chrX-46501045 | 141 | 33 | 23 | 22 | 21 |
| ZNF764-chr16-30557995 | 141 | 52 | 18 | 16 | 13 |

## Reading the machine report

`matrix.json.gz` includes every source disposition, assembly/header/program/read-group evidence, acquisition receipts, all ranked amino-acid sequences, mutation intervals, support units, independent frame/stop/sequence checks, default and primary-only results, exact-CIGAR quality sensitivity, and mitochondrial caveats.

A validated protein is a local RNA-derived window, not a reconstructed full-length protein or vaccine eligibility claim. `different_top_protein` means correct translation of RNA that differs from the reference-plus-nominated-edit expectation; inspect all checks and adjacent alleles. `protein_validation_error` is a distinct failure.

The original RNA BAMs are checksum-verified, full expanded-locus extractions. The offline corpus is deliberately selected, and its counts must not replace this matrix. Full supporting-name sets are retained losslessly via each row's `supporting_read_name_table` and each protein's `supporting_read_name_indices`, with subset SHA256. `report.protein_supporting_read_names(row, protein)` decodes and verifies them.
