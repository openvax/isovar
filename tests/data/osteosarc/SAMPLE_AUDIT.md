# Isovar sample-by-sample RNA and protein audit

## Bottom line

**All 12 sample/locus combinations now finish without crashes. All nine
with alternate reads produce a top-ranked peptide matching an independently
edited reference transcript.** MAP2 has no exact RNA deletion support in
either source; bulk PIP5K1A also has no alternate reads. These are successful
zero-alt results, not translation failures.

This PR fixes [#217](https://github.com/openvax/isovar/issues/217) (spliced
deletion-normalization crash) and
[#215](https://github.com/openvax/isovar/issues/215) (splice skip misclassified
as deletion). The results below use the proposed 1.7.12 changes, **not the
still-released 1.7.11 behavior**. Before this fix the three ONT deletion
loci PIP5K1A, H1-2 and GTF3C5 aborted in read collection. See this file's Git
history for the pre-fix audit. This PR remains unmerged for user review.

## Samples and count definitions

- **Bulk T0:** BostonGene BG003082, December 2022 tumor RNA, STAR 2.7.11b.
- **ONT T1:** UCSF June 2024 tumor single-cell long-read RNA, minimap2
  2.24-r1122, tagged/deduplicated upstream.

These differ in collection time, library preparation and processing. They
are not matched technical replicates. Differences cannot be assigned solely
to platform or tumor evolution. [README.md](README.md) documents provenance.

The report uses **full downloaded regions**: 6,239 bulk and 8,615 ONT
alignment records, not the 724-read quality-enriched CI subset. The main
table excludes secondary, supplementary and duplicate-flagged records,
then uses Isovar's minimum MAPQ of 1. No base-quality threshold, weighting,
assembly or ranking default has changed.

Counts are alignment records, not independent molecules. Paired mates can
share a name; read-name counts are not cell/UMI deduplication. “Other” is an
observed allele different from both supplied alleles, not every uninterpretable
alignment. Errors remain explicit stage/error results with null counts;
the harness never skips a crashing read to manufacture a successful count.

### Primary-only Isovar ref / alt / other

| Variant | Bulk T0 | ONT T1 |
| --- | ---: | ---: |
| PIP5K1A p.Gly474fs, 1-bp deletion | 116 / 0 / 0 | 36 / 3 / 0 |
| MAP2 p.Leu867fs, 22-bp deletion | 6 / 0 / 2 | 2 / 0 / 0 |
| H1-2 p.Ala197_Lys201del, 15-bp deletion | 328 / 30 / 46 | 3174 / 738 / 134 |
| EXOC4 p.Ser34Ile | 225 / 64 / 0 | 893 / 130 / 3 |
| GTF3C5 p.Glu503_Glu506del, 12-bp deletion | 35 / 10 / 6 | 302 / 117 / 20 |
| DYNC1H1 p.Val314Ile | 237 / 177 / 1 | 1945 / 886 / 2 |

### Unmodified read-inclusion defaults on the original regional BAMs

Defaults allow secondary alignments, exclude duplicate-flagged reads and do
not exclude supplementary alignments. This is an inclusion-policy
difference, not biological change.

| Variant | Bulk T0 | ONT T1 |
| --- | ---: | ---: |
| PIP5K1A deletion | 130 / 0 / 0 | 36 / 3 / 0 |
| MAP2 deletion | 6 / 0 / 2 | 2 / 0 / 0 |
| H1-2 deletion | 332 / 35 / 65 | 3190 / 740 / 136 |
| EXOC4 SNP | 249 / 64 / 0 | 907 / 131 / 3 |
| GTF3C5 deletion | 35 / 10 / 6 | 302 / 117 / 20 |
| DYNC1H1 SNP | 239 / 177 / 1 | 1949 / 887 / 3 |

## Protein-level validation

The protein pipeline is exercised from collected alternate reads through
RNA sequence creation, transcript matching, translation, grouping and
ranking. It uses the default requested protein length of 20, coverage 2,
minimum reference prefix 10, maximum prefix mismatches 2, no overlap
assembly and top-one output. The protein creator requests **63 cDNA bases**;
the separate `VariantSequenceCreator()` diagnostic requests 90. Their
candidate counts are not interchangeable.

### Independent expectations and annotation scope

One transcript per locus is pinned from **GRCh38 / Ensembl release 87**:

| Gene | Transcript | Protein |
| --- | --- | --- |
| PIP5K1A | ENST00000368888.8 | ENSP00000357883.4 |
| MAP2 | ENST00000360351.8 | ENSP00000353508.4 |
| H1-2 (HIST1H1C in this release) | ENST00000343677.3 | ENSP00000339566.2 |
| EXOC4 | ENST00000253861.4 | ENSP00000253861.4 |
| GTF3C5 | ENST00000372108.9 | ENSP00000361180.5 |
| DYNC1H1 | ENST00000360184.8 | ENSP00000348965.4 |

Original GTF/cDNA/peptide records, source URLs and SHA256 checksums are in
[protein_reference/](protein_reference/). No annotation download is needed
in CI. These transcript versions reproduce the reported protein changes;
this is explicitly a selected-transcript test, not an exhaustive isoform
audit or a claim that release 87 is the latest annotation.

As an additional manual check, all 12 primary-only full-region cases were
run against the complete cached Ensembl 87 annotation **without a transcript
whitelist**. The top amino-acid sequences (and three empty results) were
identical to the pinned-subset results below. This does not validate every
lower-ranked isoform; the reproducible offline oracle remains transcript-specific.

The oracle independently maps genomic positions through the raw GTF exon
intervals, verifies the reference allele against cDNA, applies the exact
allele and translates with [NCBI genetic code table 1](https://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi#SG1).
H1-2 requires reverse complementation. Unedited cDNAs must first reproduce
the original Ensembl peptide sequences exactly. Frameshift translation
continues through available downstream transcript sequence until the first
stop, rather than stopping at the original CDS boundary.

Regression checks establish the exact S34I/V314I substitutions, H1-2
AAKPK deletion, GTF3C5 EEEE deletion and PIP5K1A/MAP2 frameshifts. The
predicted full mutant reference proteins have lengths 494/906 for the two
frameshifts and 208/522 for H1-2/GTF3C5. Those full-length sequences are
**reference-derived expectations, not full-length RNA reconstructions**.
Equivalent deletion placements in GTF3C5's glutamate repeat yield the same
protein; the genomic junction is not necessarily the rightmost HGVS
protein-label junction.

For each actual RNA translation we independently check the reading frame,
amino acids, mutant interval, frameshift status and stop handling. A second,
separate comparison asks whether it equals the reference transcript with
only the nominated mutation. Expected amino acids never come from Isovar
or Varcode consequence prediction.

### Top peptides from primary-only full-region reads

Each returned sequence below matches its independently expected window.
Numbers in parentheses are supporting **read names**, not molecule counts.

| Variant | Bulk T0 top peptide | ONT T1 top peptide |
| --- | --- | --- |
| PIP5K1A | No alt reads | `RSGSSFSRRAAPVATPALLT` (2) |
| MAP2 | No alt reads | No alt reads |
| H1-2 | `KSAAKAVKPKVVKPKKAAPK` (21) | `KSAAKAVKPKVVKPKKAAPK` (628) |
| EXOC4 | `TLSTIDDVEDRENEK` (34) | `LISVIRTLSTIDDVEDRENE` (102) |
| GTF3C5 | `YESGEDEEDEEEEEDFKPSD` (7) | `ESGEDEEDEEEEEDFKPSD` (94) |
| DYNC1H1 | `KRFHATISF` (121) | `LKHGKRFHATISFDTDTGLK` (766) |

All nine supported combinations also return expected top peptides with
unfiltered read-inclusion defaults and in the coverage-1 diagnostic.
Coverage 1 can change the selected window: bulk DYNC1H1 becomes
`LKHGKRFHATISF` (13 aa, 120 names), and ONT GTF3C5 gains its leading Y.

### What this does NOT establish

- **Alt support does not guarantee the requested peptide length.** Bulk
  EXOC4's top window is 15 aa; bulk DYNC1H1's is only 9 aa; ONT GTF3C5's
  is 19 aa. None ends at a stop codon. The available compatible RNA context
  and support-based ranking limit these windows. A downstream vaccine
  workflow needing longer flanks must check length explicitly.
  **For bulk DYNC1H1, longer evidence actually exists:** requesting all
  ranked proteins yields the expected 20-aa `LKHGKRFHATISFDTDTGLK` with
  116 supporting names, but it ranks 14th behind the 9-aa window's 121
  names. Both have zero reference mismatches. Support ranks before length,
  and default top-one output hides the longer result. Even requesting 40 aa
  leaves the same 9-aa window first. This is not absence of longer RNA
  evidence; it is a concrete instance of the short-output problem tracked
  in [#90](https://github.com/openvax/isovar/issues/90). Ranking is unchanged
  in this splice-fix PR.
- **Not every translation matches the single-edit expectation.** At
  default coverage the ONT H1-2 translations have 4/15 matching windows,
  GTF3C5 3/5 and DYNC1H1 15/26; bulk DYNC1H1 has 28/29. Their top-ranked
  proteins match, and all RNA translations pass the independent frame/
  translation checks. Additional read-derived differences are not
  automatically sequencing errors or validated extra variants.
- **PIP5K1A remains weakly supported:** three ONT alternate alignments,
  two read names supporting the selected protein, no bulk alternate support.
  Successful reconstruction is not a confidence estimate.
- This does not validate every isoform, full-length expressed proteins,
  final Isovar result filters, vaxrank ranking, vaccine construct sequences,
  immunogenicity or clinical efficacy.

## Fixed bugs and remaining interpretation limits

### Spliced indels: #217 and #215

The old normalizer treated missing aligned-pair query coordinates as
deletions, conflating CIGAR D with N. MD-tagged ONT introns contain no
reference base for that path's `.upper()`, causing the run-aborting crash.
Bulk lacks MD tags despite its filename and bypassed the failing path.

Normalization now walks explicit CIGAR events and the compact
MD-reconstructed exonic/deletion reference, never expanding long introns.
Only I/D can generate indels; clipping is not insertion evidence. Shifting
stops at introns and other non-aligned events. Separately, a queried allele
overlapping a reference skip is rejected, including insertions at intron
boundaries; unrelated introns do not invalidate real exon evidence.
These distinctions follow the
[SAM specification](https://samtools.github.io/hts-specs/SAMv1.pdf).

The six former expected failures are now passing regressions. Synthetic
boundary tests cover both strands, missing/present MD, true indels on either
side of introns, adjacent D/N events, =/X operations, soft clipping, query
indices and a 10-million-base skipped intron.

### MAP2: no aligned exact-alt RNA support, not proven absent expression

The upstream export reports 135 alternate T0 tumor-WES observations.
Neither RNA source has the exact aligned deletion; coverage is sparse and
many bulk reads are soft clipped. Alignment/representation and sample
differences warrant investigation before concluding absent expression.
Nearby distinct MAP2 variants are not silently combined with this allele.

### Independent count discrepancies have explicit anchor-policy causes

The independent CIGAR oracle requires both outside anchors for every indel
observation. Isovar can accept the full reference/other allele without both
outside anchors. Alternate counts agree at all six loci in both sources.

- Bulk PIP5K1A has one extra ref observation (116 versus 115), as does MAP2
  (6 versus 5): a complete reference allele ends at the read boundary.
  Witness names are `A01231:539:HJ3NCDMXY:2:1216:11858:15922` (146M5S)
  and `A01231:539:HJ3NCDMXY:1:2152:26765:29215` (3S146M).
- ONT GTF3C5 has ten extra ref observations (302 versus 292) and one extra
  other (20 versus 19), missing an outside anchor. Example ref:
  `94d79222-1771-4ff5-8135-25180eb715b4_0`, ending at reference position
  133057905 exclusive. Other: `bd6d65db-489d-4499-94fa-f74c9e7af00c_0`,
  starting at 133057893 (0-based), after the upstream anchor.
- ONT H1-2 has two extra other observations (134 versus 132) whose deletions
  remove the downstream outside anchor:
  `c7dd0f79-4cac-4ba7-a023-443e3321723f_0` and
  `2fa6f714-f76c-4405-a4fb-496bfc2d3f9a_0`.

These are counting-contract differences, not additional mutant support.

## Quality and sequence-coverage sensitivity

No quality confidence model is inferred. A hypothetical Q20 gate on
independent primary ONT exact-alt SNP observations retains EXOC4 121/130
and DYNC1H1 855/886; bulk retains 63/64 and 174/177. STAR MAPQ 255 is not
an ordinary Phred-255 posterior. Indel anchor quality is not a quality for
deleted bases. MAPQ, raw qualities and cell/UMI metadata remain separate.

The independent 90-base RNA-candidate diagnostic represents 93/130 primary
ONT EXOC4 alt read names and 788/886 DYNC1H1 names at coverage 2. Coverage 1
recovers all 130/886. This is sequence/coverage filtering, not a collection
crash or a BQ cutoff. Those numbers must not be substituted for the distinct
protein-pipeline counts and shorter context shown above.

## Reproduce and review

```sh
python tests/data/osteosarc/fetch_sources.py /tmp/osteosarc-audit-source
python tests/data/osteosarc/audit.py /tmp/osteosarc-audit-source /tmp/new-sample-audit.json
python -m pytest tests/test_osteosarc_audit.py tests/test_osteosarc_proteins.py tests/test_spliced_indels.py -q
```

Acquisition and output paths must be new. Existing regional SAMs can be
reused. Original regional SAM and pinned reference checksums are verified
before analysis. The [JSON](sample_audit.json) records settings, source
hashes, all counts, RNA/translation stage counts, exact top peptides,
read-name support and independent comparisons. Stage exceptions preserve
previous successful results. Software metadata can differ between
environments; compare scientific rows as well as versions.

CI uses only the small original-read subsets, not network downloads or
pretend full-region counts from a selected fixture. The full-region report
is independently regenerated outside CI from checksum-verified originals.

Expansion remains [#218](https://github.com/openvax/isovar/issues/218):
all vaccine-included variants across all available RNA sources, now with
**required protein-level validation**, pinned transcript/assembly versions,
exact mutant-window/stop checks and explicit no-protein reasons. The
44-entry vaccine cohort (37 SNVs, seven deletions) and initial 50-product
RNA inventory need sample/assembly reconciliation and separate real
insertion/complex stress cases. This pilot does not complete that expansion.
