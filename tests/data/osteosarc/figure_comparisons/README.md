# Same-timepoint long/short RNA example: PIP5K1A

Two **complete indexed locus queries**, acquired September 16, 2026, from the
public CC0 osteosarc collection. The query is GRCh38
`chr1:151242177-151242180`, surrounding the source `151242178 AG>A` allele.
Original reads, qualities, CIGARs and tags are unchanged. No downsampling,
alignment editing, realignment or allele-based selection. Primary derivatives
exclude flags 256, 1024 and 2048. Full source URLs, acquisition receipts,
BAM/index SHA256 and original/retained SAM-record hashes are in
`corpus/manifest.json` (these public record hashes are not sequence validation).

| Product | Original regional alignments | Primary Isovar ref/alt/other objects | Protein, assembly on/off |
| --- | ---: | --- | --- |
| UCSF T1 tagged ONT | 95 | 55 / 4 / 1 | 46 / 46 aa |
| UCSF T1 Illumina | 182 | 129 / 1 / 0 | none / none |

The ONT protein passes independent frame, translation and reference-plus-edit
checks using the existing offline Ensembl 87 corpus models. No retained
non-focal CIGAR indels were found in the selected protein's reconstructions
(alignment-aware frame transfer shipped in 1.16.0). The Illumina product has one alternate object,
below the unchanged two-read sequence-coverage floor. Its absence of output
is **insufficient alternate support**, not evidence that the variant is absent.

Both are catalogued T1, June 2024, but different processed libraries. These are
not matched technical replicates or proven independent molecules. Cell
composition, depth, alignment and quality affect the comparison; it does not
isolate read length or prove that short reads cannot reconstruct this allele.
Do not pool their counts or infer a malignant-cell assignment from barcodes.
Source definitions: [osteosarc data catalogue](https://osteosarc.com/data/),
[variant](https://osteosarc.com/variant/PIP5K1A-chr1-151242178/).

## Reproduce acquisition

From the repository root, using HTTPS-enabled samtools and Python with pysam:

```python
import gzip, json, tempfile
from pathlib import Path
from tests.data.osteosarc.expansion.acquire import acquire_regions
from tests.data.osteosarc.figure_comparisons.rebuild import rebuild

matrix = json.loads(gzip.decompress(Path(
    "tests/data/osteosarc/expansion/audit/matrix.json.gz").read_bytes()))
variants = [v for v in matrix["variants"] if v["gene"] == "PIP5K1A"]
acquisition = tempfile.mkdtemp(prefix="isovar-pip5k1a-acquisition-")
for source in matrix["sources"]:
    if source["source_id"] in {"53f498a544883d51", "c89442609fffd3f1"}:
        receipt = acquire_regions(source, variants, acquisition, assembly="GRCh38")
        assert receipt["status"] == "ok"
rebuild(acquisition, "/tmp/pip5k1a-new-fixtures")  # destination must not exist
```

The acquisition tool downloads indexes and bounded regions, never the whole
remote BAM. Compare retained record hashes and header/allele identity as well
as archive checksums (BAM compression can vary by tool version).
Tests and figure generation use only the pinned offline derivatives.

## Context, phasing and haplotype extension

`corpus/context-evidence.json.gz` is separately pinned by
`context-manifest.json`. It retains original SAM records for both ZNF436 RNA
products, both CD109 RNA products, and five MAP2 DNA/RNA products, along with
the 35-product MAP2 audit summary and source hashes. Counts refer to products,
not independent cohorts; tagged/deduplicated ONT copies are never added.

- ZNF436: three alternate CB/UMI groups per T1 library, 3 ONT versus 8 Illumina
  read objects; 49 versus 30 aa with assembly. Even the optimistic union of
  short-read flanks is only 113 nt, below the 147 nt needed for 49 aa. The
  matched timepoint does not control library/cell composition.
- CD109: 21 ONT templates carry both alternate alleles, 207 spliced nucleotides
  apart. The short-read product has 22 focal-alt and 19 linked-alt templates,
  but no read/template or CB/UMI observes both. Three original ONT reads exactly
  support the displayed 222-nt two-edit interval, Q20 at every base. Its local
  translation is conditional on the annotated frame, not a whole-protein
  Isovar reconstruction; one 25-mer cannot include both changed codons.
- MAP2: the 22-nt listed deletion, 28-nt deletion alone, and 28-nt deletion
  plus C>A/T>G substitutions are separate explicit hypotheses. Exact anchored
  compound sequence has 24/16 T1/T2 DNA, 4/2 bulk RNA and one deduplicated T2
  ONT template at Q20. No exact isolated 22-nt sequence appears in the audited
  35 products. This is not proof of absence at arbitrary sensitivity.

`build_context.py` packages the read-only diagnostic audit and original
records; it requires `--diagnostics`, `--reference-models` (the original
Ensembl 87 model archive), and `--output`. `verify_context()` independently
recounts the pinned MAP2/CD109 observations, checks the ZNF436 molecule groups
and translations, and checks the 222-nt CD109 witnesses before plotting.
The remaining 30 MAP2 products retain audit summaries/digests, not raw SAMs.

Generate the entire gallery with
`python -m examples.osteosarc_context_figures --output-dir figures/osteosarc`
(install `isovar[plot]`). Original assembly examples come first, followed by
context/haplotype comparisons and explicit fusion RNA windows. The output is
UTC-stamped, with individual white 600-dpi PNG/SVG panels, per-example vector
PDFs, evidence JSON, a page index, and `isovar-all-figures.pdf`.

## NR2F2: matched DNA/RNA audit, September 16, 2026

The selected CeGaT RNA example contains the focal GRCh37 chr15:96875528 G>T
and a downstream deletion of GTG at chr15:96875577-96875579. A complete indexed
query of chr15:96875200-96875850 was acquired separately from the original
CeGaT blood-DNA (P116686_1), tumor-DNA (P116686_2), and tumor-RNA (P116686_3)
products. The [pinned sample metadata](https://gitlab.com/slowkow/osteosarc.com/-/blob/deaf7290a5dfa9d8d7c8ab9da5d001f69fcc2d47/scripts/data/bam-metadata-consolidated.tsv)
labels these T0 material collected December 16, 2022,
sequenced March 22, 2023. DNA headers specify hg19 (normal chrY PAR masking is
irrelevant here); RNA uses hg19. No coordinate liftover was used for counting.

`nr2f2.py` classifies the exact sequence between eight-base flanks, independently
of CIGAR placement within a repeat. The entire hg19 genomic window from the
[UCSC sequence API](https://api.genome.ucsc.edu/getData/sequence?genome=hg19;chrom=chr15;start=96875510;end=96875600)
agrees with the pinned Ensembl ENST00000394166 cDNA. Base quality must be >=20
throughout the anchored window, MAPQ >=20 (255 retained, not interpreted as a
calibrated confidence), with unmapped/secondary/supplementary/QC-fail/duplicate
flags excluded. Mates are counted once per (RG, QNAME); conflicts and other
sequences are not relabeled reference. These are template counts, not UMIs.

| Sample | GTG retained | GTG deleted | Other |
| --- | ---: | ---: | ---: |
| Blood DNA | 137 | 0 | 2 |
| Tumor DNA | 323 | 0 | 3 |
| Tumor RNA | 23 | 5 | 0 |

Four RNA templates directly phase focal T with the deletion on a single
Q20-passing segment. Tumor DNA instead has 16 directly phased T + retained-GTG
templates. All five RNA deletion templates share `81M3D17M`, the same start,
strand, end and mate start: **one endpoint family**, not five established
independent molecules. At Q30 only three deletion templates remain; lowering
MAPQ to 1 does not change deletion counts. This is **RNA-supported,
DNA-unconfirmed** context. PCR/alignment artifacts and RNA-level effects are
not distinguished; no germline or somatic origin is assigned, and zero DNA
support is not proof of absence at arbitrary sensitivity.

`corpus/nr2f2-evidence.json.gz` retains every original regional SAM record,
source headers, URLs, regional-BAM hashes, genomic reference response, and
counts at all three thresholds; `nr2f2-manifest.json` pins its SHA256.
For explicit reacquisition, fetch each `SOURCES` BAM index (`.bai`, not
`.bam.bai`), use `samtools view --no-PG -b -M -X URL INDEX
chr15:96875200-96875850 -o SAMPLE.bam`, and save the UCSC response as
`hg19-reference.json`, and the pinned `METADATA_URL` table as
`source-metadata.tsv`. Then run:

```sh
python -m tests.data.osteosarc.figure_comparisons.nr2f2 \
  --inputs /path/to/new-regions --output /path/to/new/nr2f2-evidence.json.gz
```

The destination must not exist. Tests recount original records offline; the
gallery adds separate DNA/RNA and direct-haplotype panels with caveats in the
side margin. This audit is distinct from the three-read selected translation
fixture and does not relabel that fixture as a sample-abundance estimate.

### Independent RNA-product follow-up

`nr2f2-libraries.json.gz` adds complete bounded NR2F2 queries from 23 tumor RNA
products (including the original CeGaT product), with two additional BostonGene
vendor BAMs explicitly **not assessed** because the pinned inventory has no
index. BostonGene reprocessed products are assessed. These are **25 products,
not 25 independent libraries**, and counts must not be added across rows.
The same pinned metadata table selects tumor RNA, plus the T1 PacBio mapping,
T1 Illumina source and three tagged ONT predecessors of the deduplicated ONT
products. Pooled blood data with unresolved donor attribution are not included.

UCSC chain 16 maps the entire GRCh38 chr15:96332281-96332371 interbase window
to GRCh37 chr15:96875510-96875600 in one forward, ungapped block. The independently
fetched hg38 sequence equals the original hg19 sequence across all 90 bases.
The deletion is GRCh38 chr15:96332348-96332350 (GTG), not the focal SNV.
Queries cover chr15:96331971-96332621 or chr15:96875200-96875850 as appropriate,
after verifying the BAM reference assembly. All original regional SAM records,
headers, index receipts, reference responses and sample metadata are retained.

Representative products below use the same eight-base anchors, Q20, MAPQ20,
flag exclusions and (RG, QNAME) template counting as the original audit.
"Other" includes nonmatching anchored sequences rather than relabeling them
reference. Zero callable templates means unassessable, not negative.

| Product / processing family | GTG retained | GTG deleted | Other |
| --- | ---: | ---: | ---: |
| T0 CeGaT | 23 | 5 | 0 |
| T0 BostonGene BG003082, reprocessed | 18 | 0 | 0 |
| T0 Personalis, oncoanalyser | 17 | 0 | 1 |
| T1 BostonGene BG009368, reprocessed | 85 | 0 | 0 |
| T1 Tempus, vendor | 6 | 0 | 0 |
| T2 UCLA SARC0277, reprocessed | 5 | 0 | 0 |
| T1 UCSF Illumina, CellRanger | 24 | 0 | 0 |
| T2 UCSF Illumina, CellRanger | 30 | 0 | 0 |
| T3 UCSF Illumina, CellRanger | 2 | 0 | 0 |
| T3 UCSF CD45-negative Illumina | 46 | 0 | 0 |
| T1 UCSF ONT, deduplicated | 8 | 0 | 0 |
| T2 UCSF ONT, deduplicated | 22 | 0 | 0 |
| T3 UCSF ONT, deduplicated | 0 | 0 | 0 |
| T1 PacBio | 0 | 0 | 0 |

No other acquired product supports the exact deletion at Q20, Q30, Q10, or
Q20/MAPQ1. An explicitly **sequence-only** diagnostic, including missing base
qualities, also finds no deletion outside CeGaT. PacBio has 19 reference-window
templates in that diagnostic, but all 45 regional records lack QUAL, so none
qualifies as Q20 evidence. T3 ONT has one reference-window template at Q10,
none at Q20. These outcomes do not establish absence at arbitrary sensitivity.
No other product has a Q20 focal-T call either; different sampled material,
allelic expression, assay sensitivity and processing remain relevant.

The original CeGaT five Q20 deletion templates still form **one endpoint family**
(seven templates at Q10). The deletion therefore remains **not independently
confirmed in another RNA library**, RNA-supported and DNA-unconfirmed. This
does not identify a germline/somatic origin, distinguish PCR/alignment artifacts
from RNA-level effects, or establish independent molecules.

Processing relationships are checked rather than inferred from display labels:
the two Tempus vendor deliveries share all 146 regional read names; the T1 UCSF
Illumina/CellRanger products share all 805; T1 BG009368 and the product labelled
T1 UCLA oncoanalyser share 402 of 405 each, with BG009368 in the latter's read
group. UCLA T2 products similarly overlap, and ONT deduplicated products are
subsets of their tagged predecessors. Original metadata explicitly calls the
Tempus deliveries duplicates. Different sequencing platforms or renamed read
groups are not proof of different input molecules; barcode families are not
pooled across products. CeGaT, T0 BostonGene and T0 Personalis have different
instrument/run/flowcell prefixes, supporting distinct sequencing runs, not
necessarily distinct tissue extractions or independent biological replicates.

Reproduce with the pinned `METADATA_URL` table saved as `source-metadata.tsv`:

```sh
python -m tests.data.osteosarc.figure_comparisons.nr2f2_libraries \
  --output /path/to/new/nr2f2-libraries \
  --source-metadata /path/to/source-metadata.tsv
```

Optional `--index-cache` reuses only URL- and checksum-verified original indexes
under `*/inputs/alignments/*/source-index.bai`. No whole remote BAM is fetched.
The manifest pins the resulting archive, and offline tests independently
recount all five policies from its original records. Neither acquisition nor
external data access runs in tests.

## Read-identity correction (Isovar 1.17.0, #264)

Across the 49-case expanded corpus, top protein sequences and outcome labels
are unchanged. Seven default-mode cases change counts or intermediate stage
counts (ATRX, GLIS3, NAV2, PIP5K1A, RNF213, SLC25A12, ZNF674); all primary-only
results are unchanged. NAV2's top protein support drops from 23 alignment-counted
reads to 18 sequenced segments. Two high-depth KTN1 negative controls also lose
duplicate alignment counts. No original read, reference or independent
translation expectation was modified.

The older six-locus corpus exposes a separate, genuine output change for
`bulk_star_t0/H1-2`: read2 of template `...2382:6460:24392` has competing
`99M15D46M4S` (primary) and `118M31S` (secondary) placements with different allele
calls. These observations are now uncertain rather than definitive deletion
support. One unambiguous overlapping pair (two segments, one template) remains,
so the unchanged default two-observation coverage floor yields no protein.
An explicit coverage-one diagnostic still independently validates translation;
it is not the default, and this selected fixture is not a whole-sample verdict.
