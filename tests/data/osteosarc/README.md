# Real tumor RNA at osteosarc.com cancer-mutation loci

724 original alignments at six reported somatic mutation loci: 482 from
BostonGene bulk tumor RNA aligned with STAR, and 242 from UCSF nanopore
single-cell long-read tumor RNA aligned with minimap2. These complement the
germline/alignment-edge-case examples in `../real_rna`; no reads, qualities,
CIGARs, names, headers or tags were synthesized or edited. No production
base-quality, weighting, translation or clinical policy is changed here.

## Sources and reuse

Sid Sijbrandij's osteosarcoma dataset was accessed on 2026-09-07 from
https://registry.opendata.aws/sid-osteosarc/.

The official AWS Open Data Registry specifies **CC0 1.0**. Its original entry
is retained as `source_registry.yaml`. The data were intentionally made
public; only mutation-centered RNA alignments and relevant source metadata
are included, not unrelated clinical records.

- Dataset and access instructions: https://osteosarc.com/data/
- Indexed alignment catalogue: https://osteosarc.com/bams/bams.json
- Source mutation counts: https://osteosarc.com/variants/variant_vafs_long.tsv
- Column definitions: https://osteosarc.com/variants/variant_vafs_long.columns.tsv
- License snapshot: https://raw.githubusercontent.com/awslabs/open-data-registry/main/datasets/sid-osteosarc.yaml

`selection.json` pins complete source BAM URLs and exact original GRCh38
VCF-style alleles. `source_variant_counts.tsv` retains all upstream count
rows for these six variants, including DNA tumor/normal evidence. Those
counts and protein labels are **source reports**, not our independent
RNA oracle, a certified somatic truth set, or validated peptide sequences.
Vaccine inclusion and immunogenicity claims are not used as test truth.

## Mutation examples

Counts below are exact, directly CIGAR-aligned alternate observations among
primary, non-duplicate, non-supplementary alignments **before downsampling**.
They count alignments, not independent fragments or UMI-collapsed molecules.
Coordinate anchors remain in the listed alleles; extraction removes the
shared anchor for deletions. Protein labels are copied from the source.

| Gene / source protein label | GRCh38 allele | Bulk T0 alt observations | ONT T1 alt observations |
| --- | --- | ---: | ---: |
| PIP5K1A p.Gly474fs | chr1:151242178 AG>A | 0 | 3 |
| MAP2 p.Leu867fs | chr2:209694768 CCTGGGCTACTGTGTGTTCAATA>C | 0 | 0 |
| H1-2 p.Ala197_Lys201del | chr6:26055824 ACCTTGGGCTTAGCGG>A | 30 | 738 |
| EXOC4 p.Ser34Ile | chr7:133274996 G>T | 64 | 130 |
| GTF3C5 p.Glu503_Glu506del | chr9:133057893 GGAGGAGGAGGAA>G | 10 | 117 |
| DYNC1H1 p.Val314Ile | chr14:101980529 G>A | 177 | 886 |

Variant pages are under `https://osteosarc.com/variant/{variant_id}/`;
H1-2 uses `H1_2-chr6-26055824` as its identifier. Frameshift labels alone do
not specify a complete altered protein or establish its expression.

The MAP2 22-base deletion has source-reported T0 tumor-WES support (135
alternate observations), but no exactly aligned deletion in these two RNA
regions. Only seven bulk and two ONT primary observations have both required
anchors; many bulk reads are soft clipped. This is a useful **no direct RNA
support** example, not proof of absent mutant expression. Local realignment,
alternative indel representations, sparse coverage and sample/timepoint
differences require further investigation. Nearby distinct MAP2 calls must
not be conflated with this allele. Likewise, PIP5K1A's T0/T1 contrast is not
by itself proof of biological acquisition.

## Alignment and quality caveats

- `bulk_star_t0`: BG003082, catalogued as T0 tumor RNA, December 2022,
  reprocessed with STAR **2.7.11b**, GRCh38 annotation v47. The source header
  specifies STAR and samtools 1.19.2; its read group has ID/SM but no PL or
  instrument model, so an exact sequencing instrument is not asserted.
  Despite `.md.bam` in its filename, these reads **lack MD tags**. NH/HI and
  primary/secondary alignments are preserved. MAPQ 255 corresponds to NH=1
  here; it is not a calibrated Phred-255 confidence.
- `ont_t1`: T1 UCSF tumor ONT, June 2024, **single-cell long-read RNA**, not
  the old NA12878 direct-RNA library. Header records minimap2 **2.24-r1122**,
  `-ax splice -uf --MD`, samtools 1.21, upstream secondary removal and a
  merge; source file is tagged/deduplicated. MD, CB/UB and raw barcode/UMI
  tags remain intact. Exact pore, chemistry and basecaller are not inferred
  from these alignment headers. Do not count it as raw independent molecule
  evidence or pool it with bulk T0 as a matched technical replicate.

Before downsampling, Q20 retains **855/886** exact alternate DYNC1H1 ONT
observations here, versus **4/106** alternate-base observations at the older
NA12878 ONT SNP in `../real_rna`. The latter is a base-only rather than pure
allele statistic. These different loci, libraries, processing histories and
times cannot calibrate either platform, but show why neither a universal
Q20 policy nor one fixed assumption for all ONT data is established by a
single example. Low-quality cancer-allele observations with high MAPQ are
explicitly retained for SNPs in both new sources.

For deletions, the quality audit uses the two observed flanking bases;
deleted bases have **no query quality**. Anchor quality is not deletion
confidence. `source_primary_alt_quality_retention` is descriptive only; it
does not implement or endorse a filtering or weighting formula. Site count
filters need not match this independent CIGAR walk.

## Selection and independent expectations

We fetched only indexed regions from 6.30 GiB and 34.76 GiB BAMs, yielding
6,239 bulk and 8,615 ONT regional records (about 25 MB of SAM). The fetch
interval for an anchored variant is `chrom:(pos-1)-(pos+len(ref))` in samtools
1-based coordinates. `samtools view -M` avoids repeated interval emission.

For each source/locus, choose the first **48 template names by SHA256(name)**
among records with an aligned base in that interval, without looking at
allele or quality. Then add the **two lowest-quality primary alternate-supporting
template names** at each locus, ranked by minimum query/anchor quality with
SHA256(name) tie-breaking. This explicit stress supplement catches low-BQ
cancer observations that a small blind subset missed. Retain every available
regional SAM record for the union of names; do not fetch off-region mates or
deduplicate records. Template lists for both selections are recorded in the
manifest. **The final selected collection is not a VAF estimator.**

`rebuild.py` writes reproducible gzip streams with original SAM headers and
all retained records. The independent CIGAR walker in `tests/real_rna_helpers.py`
imports no isovar code and distinguishes D from N, records exact alleles,
SAM-oriented query offsets and original quality indices, and requires both
anchors for indels. It does not repeat-normalize or perform local realignment.
The manifest pins uncompressed SAM and original regional-source checksums,
per-record hashes with multiplicity, observed alleles/qualities, and the
pre-downsampling audit. Full-source audit numbers are reproducible via the
manual recipe, not recomputed from the small CI fixtures.

## Reproduce and test

Requires curl, HTTPS-enabled samtools, and Python with pysam. The first
command downloads only indexed regions and small metadata files; it requires
a **new** destination directory. Metadata checksum drift fails explicitly
and needs review. There are no downloads during tests.

```sh
python tests/data/osteosarc/fetch_sources.py /tmp/osteosarc-new-source
python tests/data/osteosarc/rebuild.py /tmp/osteosarc-new-source /tmp/osteosarc-rebuilt
python -m pytest tests/test_osteosarc_rna.py -q -rx
```

Compare generated `.sam.gz`, `manifest.json`, `source_variant_counts.tsv`,
`source_variant_counts.columns.tsv`, `source_registry.yaml` and
`source_checksums.json` byte-for-byte against this directory.

The tests cover integrity, chromosome identity, exact allele/query-quality
extraction, indexed SNP partitioning at several MAPQ thresholds, low-BQ
cancer alleles, preserved aligner/cell/UMI metadata, source allele definitions,
and deletion-to-sequence creation. Real ONT deletions expose a **new crash**:
[isovar #217](https://github.com/openvax/isovar/issues/217), where normalization
calls `.upper()` on the missing reference base of a CIGAR N. Three strict
expected failures recognize only that specific failure path; unrelated
errors or bad fixtures still fail. This is separate from the false-deletion
splice-skip bug [#215](https://github.com/openvax/isovar/issues/215). Both must
be fixed before describing spliced indel handling as robust. These fixtures
are not end-to-end translation, somatic calling or clinical validation.
