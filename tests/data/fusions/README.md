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
original T1/T2 ONT deduplicated BAMs linked in each receipt. RNA event
coordinates follow the catalogue's `scripts/fusions/fusions.tsv`, not its
reference-only canonical transcript drawings. TPST1 has different DNA and RNA
breakpoints; this fixture uses the RNA junction. FOXO3 and PARD3B use the
listed DNA junction with direct RNA confirmation. PARD3B retains the observed
12-nt `TGTGCTGATATG` insertion.

`build_osteosarc.py --alignments /path/to/regional-bams --output corpus` selects
primary records with original SA tags, excludes secondary/duplicate/QC-fail
records, requires MAPQ >=20 for both partner mappings and Q>=10 across every
retained base, and retains unique physical reads/cell-UMI groups. It never
fills RNA bases from reference. It starts with 60-nt flanks, then tries 45,
30, and 20 only when no exact group has two distinct fragments. This is an
illustrative group-selection rule, not a sensitivity or discovery benchmark.

The regional BAMs were acquired with indexed `samtools view --no-PG -b -M -X`
from the original public objects; full commands and BAM/index hashes are
retained in each fixture. Regions are chr7:66205000-66206000,
chr7:66127200-66128200, chr6:108609245-108613245,
chr17:63743980-63747980, chr2:205308103-205312103, and
chr9:22005648-22009648. The builder requires the complete cached Ensembl 87
annotation/cDNA, not a subset registered under a generic genome name.

All five retained windows remain `unresolved_frame` with the supplied
annotation. TPST1 joins before its annotated CDS start; FOXO3 and PARD3B
windows have no exact collinear donor transcript match. No repeated PARD3B
T2 window passes this specific extraction rule; that is not absence of RNA
support. Tests independently compare every retained sequence and quality to
the original SAM record, and synthetic coding examples independently check
translation with the published NCBI table. No coding protein is fabricated
for these real candidates.
