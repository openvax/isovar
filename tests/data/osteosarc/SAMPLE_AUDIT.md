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
as deletion). The results below use the proposed 1.8.0 changes, **not the
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
then uses Isovar's minimum MAPQ of 1. No base-quality threshold, weighting
or assembly default has changed. Protein context generation and ranking
now use the balanced policy described below.

This audit uses `ReadCollector()` (unmerged overlapping mates); `run_isovar()`
normally enables mate merging, which can change record/candidate counts.
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
ranking. The new default targets **49 aa for 25mer design**, with coverage 2,
minimum reference prefix 10, maximum prefix mismatches 2, no overlap
assembly and top-one output. Among candidates retaining at least **90% of
the best mutant candidate's compatible read-name support**, it maximizes
the number of full mutation-containing 25mer placements. If none is
available within budget, it prefers useful partial context without quietly
lowering the threshold. Read-name support is not confidence or evidence
that every name spans the complete peptide.

A bounded RNA-context ladder (20, 25, 29, 33, 37, 41, 45, 49 aa targets)
keeps shorter reference-matching candidates available when long noisy
flanks fail. The longest requests **149 cDNA bases**. Identical candidates
are deduplicated without inventing additional reads or joining conflicting
haplotypes. The separate `VariantSequenceCreator()` diagnostic requests 90
bases; its candidate counts are not interchangeable with this pipeline.
See [policy and configuration](../../../PROTEIN_SELECTION.md).

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
| PIP5K1A | No alt reads | `VFKKIPLKPSPSKKFRSGSSFSRRAAPVATPALLTSHRSLGNTRHK` (2) |
| MAP2 | No alt reads | No alt reads |
| H1-2 | `SPKKAKVAKPKKAAKSAAKAVKPKVVKPKKAAPKKK` (21) | `KKAAKSAAKAVKPKVVKPKKAAPKKK` (569) |
| EXOC4 | `YRSTVSKSKDPSGLLISVIRTLSTIDDVEDRENEKGRLEEAYEKCDRDL` (32) | `GLLISVIRTLSTIDDVEDRENEKGR` (97) |
| GTF3C5 | `FSSSAKADGGKEQLTYESGEDEEDEEEEEDFKPSDGSENEME` (7) | `KEQLTYESGEDEEDEEEEEDFKPSDGSEN` (90) |
| DYNC1H1 | `RESPEVLLTLDILKHGKRFHATISFDTDTGLKQALETVNDYNPLMKD` (111) | `LLTLDILKHGKRFHATISFDTDTGLKQALETVN` (692) |

All nine supported combinations also return expected top peptides with
unfiltered read-inclusion defaults and in the coverage-1 diagnostic.
Coverage 1 retains these selected amino-acid sequences but can change
their supporting names (ONT H1-2 583, EXOC4 101, DYNC1H1 702).

### Context versus retained support

Each entry is **length / selected names / best candidate names / valid
25mer placements**, using primary-only reads and the balanced default.

| Variant | Bulk T0 | ONT T1 |
| --- | ---: | ---: |
| PIP5K1A | No alt | 46 / 2 / 2 / 21 |
| MAP2 | No alt | No alt |
| H1-2 | 36 / 21 / 21 / 12 | 26 / 569 / 631 / 2 |
| EXOC4 | 49 / 32 / 34 / 25 | 25 / 97 / 102 / 1 |
| GTF3C5 | 42 / 7 / 7 / 17 | 29 / 90 / 94 / 5 |
| DYNC1H1 | 47 / 111 / 121 / 23 | 33 / 692 / 767 / 9 |

49 aa is sufficient for all 25 placements **only for a centered
single-residue mutation**. A mutation-free window is never counted. Pure
deletion junctions require strict overlap across the new adjacency; wider
edits and frameshifts can need more context to cover every placement.

DYNC1H1 p.Val314Ile is well inside a 4,646-aa protein. Its old 9-aa top
result was not a protein boundary or stop: the 20-aa support-first request
and ranking preferred 121 compatible names for `KRFHATISF` over 116 for
`LKHGKRFHATISFDTDTGLK`. Increasing the target alone still ranked the short
window first. Some of the extra names disagree outside that short core;
we must not simply assign their support to a longer conflicting sequence.

| DYNC1H1 preference | Bulk length / names / 25mer placements | ONT length / names / 25mer placements |
| --- | ---: | ---: |
| Historical support-first, explicit target 20 | 9 / 121 / 0 | 20 / 766 / 0 |
| Balanced, retain 95%, target 49 | 35 / 115 / 11 | 25 / 744 / 1 |
| Balanced, retain 90%, target 49 (default) | 47 / 111 / 23 | 33 / 692 / 9 |
| Context-first, target 49 | 47 / 111 / 23 | 49 / 516 / 25 |

All those DYNC1H1 sequences match their independently expected windows.
The 49-aa ONT context-first sequence is
`EKRESPEVLLTLDILKHGKRFHATISFDTDTGLKQALETVNDYNPLMKD`.
The old 20-aa and new fallback RNA budgets differ by one base, hence ONT's
766 historical names versus the new ladder's best 767. The default bulk
47-aa window retains 91.7% of candidate support, ONT's 33 aa 90.2%; these
are not fractions of every alternate alignment or calibrated confidence.

**Context-first is intentionally not the default.** ONT H1-2 then selects
`SPKKAKVAKPKKAAKSAAKAVKPKVVKPKEGGAQEEIGERLLLKPKRLF` (49 aa, only six
names), whose tail differs from the nominated-edit expectation. The actual
RNA still passes the independent frame/translation checks; this is weak,
discordant RNA context, not proof of an additional tumor mutation. Its
reported outcome is `different_top_protein`, not silently called validated.

All five policies/diagnostics are explicit in JSON: `protein_default`,
`protein_coverage1`, `protein_support20`, `protein_balanced95`, and
`protein_context`. The relative threshold never changes allele counts or
the absolute sequence/reference filters. Vaxrank still explicitly requests
35 aa; [vaxrank #415](https://github.com/openvax/vaxrank/issues/415) tracks
coordinating its default request. Its current padding value 12 requests
49 aa for 25mers. Actual final vaccine-selection validation remains
[vaxrank #414](https://github.com/openvax/vaxrank/issues/414).

### What this does NOT establish

- **Alt support does not guarantee the requested peptide length.** All
  nine positive cases now yield at least one valid 25mer under balanced
  90%, but ONT H1-2 falls back to 20 aa under the stricter 95% policy.
  Other samples can lack any eligible full window. This addresses the
  demonstrated ranking limitation in [#90](https://github.com/openvax/isovar/issues/90),
  not every possible cause of short output. There is no reference padding.
- **Not every translation matches the single-edit expectation.** Balanced
  top-ranked proteins match in these nine cases; the JSON retains counts
  of other translations and explicit context-first disagreement. All RNA
  translations pass the independent frame/translation checks. Additional
  read-derived differences are not automatically errors or validated variants.
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
protein-pipeline counts and multi-scale contexts shown above.

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
