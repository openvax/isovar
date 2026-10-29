# Sid matched-normal reassessment

This report inventories all 357 VCF products in Osteosarc snapshot
`2f1ad541f104b680ba0ee0652603216bc40c16b1486ab837a79dd4944d640084` and assesses
208,063 nominated alleles against 15 selected normal calling/reference-confidence
products. The reviewed selection is pinned in
[sid_normal_inventory.json.gz](sid_normal_inventory.json.gz). Runtime selection
uses complete source IDs and exact VCF sample columns, not filename heuristics.

## Results

The expanded v3 assessment covers **all 251,934 previously supplied background
alleles** in both requested and acquired territory. Its [reconciliation receipt](sid_v3/reconciliation.json)
pins the native assessment SHA-256
`e563d5b3c89b71e680e0cf9529eacf9b76a2ddea063d34e64f02e938810aade2`.
Every source and exclusion is listed in the [source inventory](sid_v3/source_inventory.tsv).

| Normal state | Nominated alleles |
|---|---:|
| Confidently reference | 195,436 |
| Conflicting calls | 5,327 |
| Uncalled | 3,342 |
| Other allele observed | 2,244 |
| Nominated ALT observed | 1,714 |
| **Total** | **208,063** |

There are **7,041 origin-review flags**: 6,843 newly flagged, 198 retained and
16 previous flags no longer supported under this policy. These are normal-DNA
evidence states, not revised RNA counts or a declaration of tumor specificity.
[Every changed flag](sid_v3/origin_changes.json.gz) retains its prior/new state.

The new background contains **646,792 alleles**: 396,524 occur in at least one
of the two original normal products and 250,268 occur only in additional selected
products. Relative to the prior background, 405,430 are newly supplied and
10,572 are not retained under the new policy. The original-product audit identifies
1,114 previously supplied alleles that were listed but not genotype-selected in
either original normal product. Some can have independent support in additional
products; that audit is not itself a blanket deletion list. The companion
[new-allele](sid_v3/new_background_alleles.json.gz),
[not-retained](sid_v3/previous_background_not_retained_in_scope.json.gz) and
[uncalled ALT](sid_v3/legacy_unindexed_alternates.json.gz) files retain exact identities.

**MT-ND5, chrM:12994 G>A, remains uncalled in normal.** Two BostonGene
reference-confidence products contain no-call records with zero minimum depth;
the other selected products have no qualifying record there. Its
[full normal evidence](sid_v3/mitochondrial_evidence.json.gz) preserves those
records and reasons. This is a calling/coverage limitation, not a lost contig alias.

## Source selection

The 15 products are not 15 independent specimens. They include repeated calling
and alignment products from BostonGene, Personalis and UCLA normal inputs. The
published registry distinguishes the BostonGene January 2023 draw, UCLA June 6,
2024 draw and Personalis June 11, 2024 draw; its older `T0_blood` label combines
different draws. UCLA paired-normal realignments share an underlying library.
BostonGene products are conservatively grouped by normal-input lineage: the
registry establishes a reused draw, but that group does not certify identical
physical library preparation. Every product remains individually identifiable.
Reviewed VCF-to-specimen assignments await upstream catalogue support in
[Osteosarc #120](https://github.com/iskandr/osteosarc/issues/120).

Selection keeps both previous Strelka variant-only normal sources, additional
Strelka reference-confidence products, a variant-only HaplotypeCaller product,
FreeBayes products and DRAGEN normal products. Exclusions retain explicit reasons
for unresolved specimen/build assignments, somatic/annotation products, redundant
representations, structural/CNV scope and missing confidence fields. A filename
or VCF header declaring a field does not establish a qualified genotype.

## Evidence and background scope

The library's default policy requires PASS, a complete genotype, quality >= 20
and depth >= 5. Strelka explicitly uses GQX and DPI for literal indels; reference
blocks require MIN_DP. FreeBayes records without the declared confidence/filter
requirements remain uncalled, without substituting site QUAL. Processing copies
do not vote. ALT, conflicting calls, confidently reference, other ALT, uncalled,
variant-only absence and unavailable evidence remain distinct.

Background acquisition covers nomination intervals with 500-base flanks **and**
the union of Ensembl 87/115 exons, pinned ORFdata coding blocks and previously
supplied background REF spans. The additional union contains 466,465 regions and
217,230,718 bases. This includes distant spliced exons that local flanks would
miss. It is declared analysis territory, not a claim of complete normal-genome
confidence.

Only accepted, genotype-indexed literal ALTs enter the background. An old report
preparation bug selected every listed ALT when any ALT index was called. For
example, `G -> GA,GAA; GT=1/1` supplies GA, not GAA. Reconciliation separately
identifies these legacy entries; other differences can reflect scope or policy
and are not automatically classified as errors.

The new region preparation also fixes
[Isovar #612](https://github.com/openvax/isovar/issues/612): ORFdata coding blocks
are 1-based inclusive. The old preparation dropped their first bases, including
valid one-base blocks. The current conversion is `[start - 1, end)`, with a
regression using native FAM107B block geometry. The catalogue itself is valid.

## Reproduce

Use the original pinned cohort inputs, including the unflagged nomination list,
the previous supplied background and the ORFdata catalogue. Ensembl 87 and 115
annotation databases must already be installed. Run from the repository root:

```bash
COHORT_INPUTS=/path/to/cohort_reanalysis_20261027

python -m reports.matched_normal.prepare_sid_background_regions \
  --ncorf-catalog "$COHORT_INPUTS/human_all_ncorfs.json.gz" \
  --previous-background "$COHORT_INPUTS/inputs/Sid.germline.exonic_contexts.json.gz" \
  --output normal-scope/Sid.background-regions.json.gz

python -m reports.matched_normal.acquire_sid_normal \
  --variants "$COHORT_INPUTS/inputs/Sid.small_variants.json.gz" \
  --additional-background-regions normal-scope/Sid.background-regions.json.gz \
  --output normal-acquisition --workers 1

# Use the content-addressed manifest path printed by acquisition.
NORMAL_MANIFEST=normal-acquisition/normal_manifest.PASTE_ID.json.gz
python -m isovar normal-evidence \
  --manifest "$NORMAL_MANIFEST" \
  --output Sid.matched_normal.json.gz

python -m reports.matched_normal.reconcile_sid_background \
  --assessment Sid.matched_normal.json.gz \
  --previous-variants "$COHORT_INPUTS/inputs/Sid.small_variants.with_origin_flags.json.gz" \
  --previous-background "$COHORT_INPUTS/inputs/Sid.germline.exonic_contexts.json.gz" \
  --output-dir normal-reconciliation
```

Acquisition sends only public source URLs and public HTTP identities. Coordinates
and genotype filtering remain local. Each full original download is hashed and
size-checked in temporary storage before local parsing, then discarded. Complete
regional derivatives and indexes are pinned; transient transport failures can be
retried without accepting incomplete files or changing source identity. Completed
acquisitions are verified before reuse. No DNA alignments are acquired.

## Interpreting the output

`reconciliation.json` records all nomination state counts, source-specific
outcomes, changed origin-review flags and background comparisons. Companion
compressed files retain every changed flag, newly supplied background allele,
previous in-scope allele not retained, legacy uncalled ALT and mitochondrial
assessment. `source_inventory.tsv` lists every source, exact sample column,
declared specimen/library group, selection reason and original URL.

The native normal receipt additionally retains every selected-source assessment,
record genotype, rejection reason, acquisition interval and file/index digest.
The earlier local-flank-only v2 receipt remains unchanged. The expanded v3 result
is a new normal/background assessment; earlier RNA and protein receipts still
describe their original supplied background.

For a new cohort run, use unflagged Sid nominations and pin the new receipt as
the task's `normal_evidence`, replacing `germline_variants`. Reuse the appropriate
nomination subset for a targeted pass. Do not attach a new normal receipt to an
old completed RNA result or carry obsolete generated origin flags forward.
See the [shared API/CLI contract](../../docs/matched-normal-evidence.md).
