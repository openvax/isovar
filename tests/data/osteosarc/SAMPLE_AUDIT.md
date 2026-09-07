# Isovar sample-by-sample RNA audit

Spec: compare full downloaded regions for the existing six loci in two RNA
sources; distinguish independent CIGAR evidence from Isovar defaults and
primary-only counts; preserve crashes as errors with null counts; report
sequence creation separately from translation. Add a reproducible report and
offline tests, with no production behavior change. This PR stays open for
user review; #217 is a subsequent fix, not part of this audit.

## Bottom line

**The tested SNPs reach RNA sequence creation in both samples. Two tested
bulk-RNA deletions do too. Three ONT deletion loci abort in read collection
despite directly aligned alternate support.** MAP2 is different: neither
source provides the exact aligned deletion, and Isovar returns zero alt
reads without crashing. None of this yet validates translated peptides.

This is the behavior of released Isovar 1.7.11, re-audited on the 1.7.12
report branch (version bump only; no production behavior changes). Exact
software versions and settings are in [sample_audit.json](sample_audit.json).

## Which samples and counts?

- **Bulk T0:** BostonGene BG003082, December 2022 tumor RNA, STAR 2.7.11b.
- **ONT T1:** UCSF June 2024 tumor single-cell long-read RNA, minimap2
  2.24-r1122, tagged/deduplicated upstream.

These differ in collection time, library preparation and processing. They
are not matched technical replicates, and differences cannot be assigned
solely to platform or tumor evolution. Source metadata and the full URLs
are documented in [README.md](README.md).

The report uses the **full downloaded regions** (6,239 bulk and 8,615 ONT
alignment records), not the 724-read quality-enriched CI fixture subset.
It applies no base-quality threshold or weighting. The main table removes
secondary, supplementary and duplicate-flagged alignments before running
Isovar, and retains its minimum MAPQ of 1. It counts alignment records,
not independent molecules; paired mates can share a name. “Other” means
an observed allele that is neither the supplied reference nor alternate,
not all reads that could not be interpreted.

### Actual Isovar ref / alt / other counts, primary-only

| Variant | Bulk T0: ref / alt / other | ONT T1: ref / alt / other |
| --- | ---: | ---: |
| PIP5K1A p.Gly474fs, 1-bp deletion | 116 / 0 / 0 | **Error**; independent CIGAR sees 3 alt |
| MAP2 p.Leu867fs, 22-bp deletion | 6 / 0 / 2 | 2 / 0 / 0 |
| H1-2 p.Ala197_Lys201del, 15-bp deletion | 328 / 30 / 46 | **Error**; independent CIGAR sees 738 alt |
| EXOC4 p.Ser34Ile | 225 / 64 / 0 | 893 / 130 / 3 |
| GTF3C5 p.Glu503_Glu506del, 12-bp deletion | 35 / 10 / 6 | **Error**; independent CIGAR sees 117 alt |
| DYNC1H1 p.Val314Ile | 237 / 177 / 1 | 1945 / 886 / 2 |

For the three errors, there is **no completed Isovar count**. The stated
alternate observations come from the independent CIGAR walker, not from
continuing Isovar after dropping crashing reads. The JSON has `counts: null`
and the exception/stage. The audit harness catches failures between loci;
the production collection call itself still aborts.

### What changes with unmodified Isovar defaults?

The JSON also runs `ReadCollector()` directly on the original regional BAM,
without the primary-only prefilter. Its defaults allow secondary alignments,
exclude duplicate-flagged reads, and do not exclude supplementary alignments.

| Variant | Bulk T0: ref / alt / other | ONT T1: ref / alt / other |
| --- | ---: | ---: |
| PIP5K1A deletion | 130 / 0 / 0 | Error |
| MAP2 deletion | 6 / 0 / 2 | 2 / 0 / 0 |
| H1-2 deletion | 332 / 35 / 65 | Error |
| EXOC4 SNP | 249 / 64 / 0 | 907 / 131 / 3 |
| GTF3C5 deletion | 35 / 10 / 6 | Error |
| DYNC1H1 SNP | 239 / 177 / 1 | 1949 / 887 / 3 |

This is an inclusion-policy difference, not evidence of a biological change.
For example, bulk DYNC1H1's 177 primary alt alignments represent **122 read
names**; EXOC4's 64 represent **35 names**. Those names are still not an
explicit cell/UMI molecule count. Likewise, additional supplementary ONT
alignments must not automatically be treated as independent observations.

## Unexpected failures versus expected limitations

### 1. Confirmed crash: spliced, MD-tagged deletion normalization (#217)

PIP5K1A, H1-2 and GTF3C5 in ONT T1 fail with
`AttributeError: 'NoneType' object has no attribute 'upper'` inside
`ReadCollector._left_aligned_indel_interval_for_variant`. A CIGAR N splice
skip supplies no reference base for `.upper()`. Even an unrelated intron
before the queried locus can trigger the crash. This is not insufficient
coverage and is not fixed by a base-quality threshold.

The bulk BAM happens to **lack MD tags**, despite `.md.bam` in its filename.
That bypasses this normalization path. Its successful deletion collection
does not establish that normalization is correct across RNA sources.
Reproduction and acceptance criteria: [#217](https://github.com/openvax/isovar/issues/217).

### 2. Separate confirmed wrong-classification path: splice skip as deletion (#215)

Low-level extraction can return a CIGAR N as an empty deletion allele.
The earlier real-RNA fixtures reproduce this with actual splice-junction
negative probes. This audit does **not** establish that the six supplied
cancer alleles contain false-positive support from that path, but fixing
the crash must not leave the wrong classification intact. SAM distinguishes
reference skips N from deletions D; they cannot be interchanged.
[Issue #215](https://github.com/openvax/isovar/issues/215),
[SAM specification](https://samtools.github.io/hts-specs/SAMv1.pdf).

### 3. MAP2: no directly aligned RNA support, not proven biological absence

The upstream count export reports 135 alternate T0 tumor-WES observations
for this allele. Our independent RNA CIGAR walk finds no exact deletion in
either source; only seven bulk and two ONT primary observations have both
required anchors. Many bulk reads are soft clipped. Sparse coverage, local
alignment/representation issues and sample differences need investigation
before concluding that mutant MAP2 is not expressed. Nearby distinct MAP2
variants are not silently combined with this 22-base deletion.

### 4. Two ref-count discrepancies are explained by anchor policy

The strict independent indel oracle requires both outside anchors even for
reference-supporting reads; Isovar accepts a complete observed reference
allele that ends at the read boundary. This gives one extra primary bulk
ref observation at PIP5K1A (116 versus 115) and MAP2 (6 versus 5).

Original witnesses are PIP5K1A read
`A01231:539:HJ3NCDMXY:2:1216:11858:15922`, flag 83, CIGAR `146M5S`, and MAP2
read `A01231:539:HJ3NCDMXY:1:2152:26765:29215`, flag 163, CIGAR `3S146M`.
Both cover the full reference allele but lack its downstream outside anchor.
These discrepancies are not automatically bugs or extra deletion support.

## Allele support is not the same as an RNA candidate or translated peptide

Using `VariantSequenceCreator()` defaults (coverage 2, preferred length 90,
no overlap assembly), primary bulk H1-2 and GTF3C5 produce **4 and 2 RNA
sequence candidates** respectively. Both SNPs produce candidates in both
samples. ONT deletions never reach this stage because collection failed.

Primary ONT SNP sensitivity to the RNA sequence coverage requirement:

| SNP | Alt read names entering | Names represented at coverage 2 | Names represented at coverage 1 |
| --- | ---: | ---: | ---: |
| EXOC4 | 130 | 93 (5 candidates) | 130 (42 candidates) |
| DYNC1H1 | 886 | 788 (26 candidates) | 886 (163 candidates) |

The missing names here reappear when only the coverage setting changes.
This is coverage/sequence-compatibility filtering, not a read-collection
crash or a BQ cutoff. It may merit sensitivity analysis, but is not by itself
an unexpected bug. Coverage 1 is a diagnostic, not a recommended new default.
Candidate counts are not isoform counts or numbers of validated proteins.

The audit stops **before reference-transcript matching, translation, final
Isovar result filters and vaxrank ranking**. No conclusion about delivered
mutant peptides, vaccine correctness or efficacy follows from this report.

## Quality implications

For the primary ONT exact-alt SNP observations, a hypothetical Q20 gate
retains EXOC4 **121/130** and DYNC1H1 **855/886**; bulk retains **63/64** and
**174/177** respectively. These counts are descriptive sensitivity, not
calibration or a variant-confidence posterior. STAR MAPQ 255 is not an
ordinary Phred-255 posterior, and indel anchor quality is not a quality for
deleted bases. Raw qualities and MAPQ/cell/UMI metadata stay separate.

## Reproduce and review

```sh
python tests/data/osteosarc/fetch_sources.py /tmp/osteosarc-audit-source
python tests/data/osteosarc/audit.py /tmp/osteosarc-audit-source /tmp/new-sample-audit.json
python -m pytest tests/test_osteosarc_audit.py -q
```

The first path and output JSON must be new. If the original regional SAMs
already exist, run only the second command. Source checksums are verified
before analysis. Both default and coverage-1 sequence results are recorded,
and sequence-stage errors preserve any successfully collected allele counts.
Runtime software metadata can differ on another environment; compare counts
and status rows as well as versions. CI uses only the small checked-in reads,
never downloads the full regions and does not pretend to rerun full-source
counts from a biased subset.

Expansion is tracked in [#218](https://github.com/openvax/isovar/issues/218):
all vaccine-included variants across all available RNA samples, with explicit
ref/alt/other/error states and separately identified insertion/complex-allele
stress cases. The [site currently lists 44 vaccine-included entries](https://osteosarc.com/variants/);
their exported alleles are 37 SNVs and seven deletions. The initial RNA BAM
inventory has 50 products, not 50 independent biological samples, and needs
reconciliation with the broader [data catalogue](https://osteosarc.com/data/).
This six-locus/two-source pilot is not that exhaustive audit.

Requested sequence: review this PR, revise it, then merge; only afterward
implement the #217 fix. No merge, deployment or production bug fix is part
of this reporting pass.
