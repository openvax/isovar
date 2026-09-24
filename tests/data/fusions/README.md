# Original RNA fusion windows

`*.input.json.gz` are compact supplied-fusion inputs with original SAM records,
source acquisition receipts and hashes, read offsets, exon/cDNA/CDS models
(Ensembl 87 / GRCh38), and the explicit breakpoint contract. Uncompress before
passing an input to `isovar fusion`. Manifest hashes pin every fixture.

These are **partial directly observed RNA windows**, not whole-transcript
assemblies or evidence of coding protein expression. A repeated exact group
is selected independently within each sample; counts do not estimate total
fusion abundance and processed products are not independent samples.

Sources: [osteosarc fusion catalogue](https://osteosarc.com/fusions/) and the
original T1/T2 ONT tagged BAMs linked in each receipt. RNA event
coordinates follow the catalogue's `scripts/fusions/fusions.tsv`, not its
reference-only canonical transcript drawings. TPST1 has different DNA and RNA
breakpoints; this fixture uses the RNA junction. FOXO3 and PARD3B use the
listed DNA junction with direct RNA confirmation. PARD3B retains the observed
12-nt `TGTGCTGATATG` insertion.

`build_osteosarc.py --alignments /path/to/regional-bams --output corpus` selects
primary records with observed supplementary partners, excludes secondary/duplicate/QC-fail
records, requires reciprocal SA path links, MAPQ >=20 for both actual partner
mappings and Q>=10 across every retained base, and retains read/CB-UMI labels
(not proven independent molecules). It never
fills RNA bases from reference. It starts with 60-nt flanks, then tries 45,
30, and 20 only when no exact group has two distinct fragments. This is an
illustrative group-selection rule, not a sensitivity or discovery benchmark.

SA CIGARs are approximate in minimap2: they are used only to validate path
identity, never to map junction bases. Actual record CIGARs, hard clipping,
strand and original sequence agreement determine the mapping. Missing
partners remain unavailable; pairs and alternatives cannot supply them.

The regional BAMs were acquired with indexed `samtools view --no-PG -b -M -X`
from the original public objects; full commands and BAM/index hashes are
retained in each fixture. Regions are chr7:66205000-66206000,
chr7:66127200-66128200, chr6:108609245-108613245,
chr17:63743980-63747980, chr2:205308103-205312103, and
chr9:22005648-22009648. The builder requires the complete cached Ensembl 87
annotation/cDNA, not a subset registered under a generic genome name.

All five retained windows remain `unresolved` with the supplied
annotation. TPST1 joins before its annotated CDS start; FOXO3 and PARD3B
windows have no exact collinear donor transcript match. No repeated PARD3B
T2 window passes this specific extraction rule; that is not absence of RNA
support. Tests independently compare every retained sequence and quality to
the original SAM record, and synthetic coding examples independently check
translation with the published NCBI table. No coding protein is fabricated
for these real candidates.

## September 17 audit of #288

All 21 previously retained original windows (TPST1 5/7, FOXO3 4/3, PARD3B
2 T1) were matched to their actual tagged-BAM partner records: **their sequences
and base mappings were unchanged**. This verifies those selected windows, not
the former extraction algorithm in general. A real Sid regression independently
demonstrates an internal insertion that an approximate SA CIGAR relocates.

Rebuilding with the tagged predecessor changes the selected-group counts:
TPST1 T1/T2 61/122 CB-UMI labels (96/132 records), FOXO3 5/5 labels (7/5
records), PARD3B T1 2 labels/records. These counts must not be compared as
independent abundance estimates with deduplicated products. TPST1/FOXO3 retain
the same windows. PARD3B now has 60-nt flanks instead of the older shorter
selection, with the same observed 12-nt insertion; this is additional tagged
input support, not evidence that the old displayed bases were incorrect.

The nine-candidate RNA extension and its limitations are documented in
`../osteosarc/figure_comparisons/RNA_FOOTPRINTS.md`.
