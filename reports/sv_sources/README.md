# SV source interpretation and RNA replay

These receipts document input interpretation for Isovar #575 and #588. They do
not establish translation, tumor specificity or independent biological events.
See [the API/CLI guide](../../docs/sv-source-validation.md) for the source profiles,
original downloads, pinned author revisions and unresolved coordinate conventions.

| Original source | Rows | Resolved adjacency | Unresolved | Invalid | PCR positive |
|---|---:|---:|---:|---:|---:|
| COLO829 Paulin Paz v2 | 54 | 46 | 7 | 1 | Not recorded |
| HCC1395 SVScope supplement | 1,118 | 157 | 181 | 780 | 32 |

`COLO829.source.json.gz` and `HCC1395.source.json.gz` retain every original row,
its line number, source checksum, explicit transformations, caller identifiers,
PCR tri-state and interpretation reasons. Each report records the producing
Isovar and Varcode versions and interpretation-module hash.

## Comparison with author files and reference

`source_comparison.json.gz` was produced by `compare_original_sources.py`.
The inputs and driver are SHA-256 pinned in that receipt. The reference is the
GRCh38 full analysis set plus decoys/HLA used for the local RNA alignment,
SHA-256 `3b103f4742abfd54938fb0333e19ad067635c8eb86f1dbf0ce44b165c4292b50`.
This check uses existing reference sequence; it does not acquire DNA reads.

All four COLO829 embedded breakend alleles agree with the amended Minda source
by CHROM, POS, REF and interpreted ALT. Three have different event IDs in the
amended file. No other event is equated solely by its ID or an approximate locus.
The zero-width inversion remains invalid, and missing insertion sequences remain
unresolved.

All 1,096 shared HCC1395 rows have identical core columns between the author
GitHub ensemble and the PCR supplement. Of those, 157 rows are entirely identical
and 939 differ only in INFO. The supplement adds 22 literal consensus alleles.
None of those 22 REF sequences matches GRCh38 at POS; five match 49 bases upstream,
and 17 have no exact REF match within 500 bases either side. This is an exact
sequence check, not a realignment or correction of the source coordinates.
The author's writer uses a normal consensus as REF and a discovery-window start
as POS. Those records require validated anchoring before ordinary allele use.

Reproduce with the original files linked in the guide:

```bash
python reports/sv_sources/compare_original_sources.py \
  --hcc tests/data/sv_sources/HCC1395.svscope.pcr.vcf \
  --author-hcc HCC1395_minda_ensemble.Filter.vcf \
  --colo tests/data/colo829_whole_run/colo829-sv-grch38.vcf \
  --minda-colo colo829_benchmark_grch38_amended_vaf.vcf \
  --reference GRCh38_full_analysis_set_plus_decoy_hla.fa \
  --output source_comparison.json.gz
```

## Original-RNA replay of interpreted COLO829 breakends

`COLO829.interpreted_breakends_rna.json.gz` records the standard `sv-audit`
results for the four newly interpreted events, each in both orientations. It
uses the original SRR8615617 selected-sites BAM and unchanged discovery budgets.
The shared manifest contains all 188 source records, including ten remaining
unassessed definitions across its two inventories. Native snapshot and completion
receipt hashes were verified before this compact receipt was produced.

SV22 and SV29 produced no retained paths. SV28 and SV45 each produced eight paths
per orientation, all remaining regional candidates. All eight views had zero
direct-junction paths and zero event-linked proteins. Selected-site RNA cannot
exclude unsampled distant sequence or establish absence in the complete library.
These are additional assessed nominations, not additional protein candidates.
The receipt preserves its actual producer versions; no old result is relabelled
as a run with a newer library.
