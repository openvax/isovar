# Collecting and counting RNA reads

Isovar collects the RNA reads that overlap each variant with a `ReadCollector`.
Create one to change which reads are used. The defaults are shown:

```python
from isovar import run_isovar, ReadCollector

read_collector = ReadCollector(
    min_mapping_quality=1,
    use_duplicate_reads=False,
    use_secondary_alignments=True,
    use_soft_clipped_bases=False,
    # Merge overlapping mates of one fragment into a single observation.
    merge_overlapping_fragments=True,
    # Keep reads without QUAL; their base qualities stay unknown.
    use_reads_without_base_qualities=True,
    # Optional predicate on each original pysam record.
    read_filter=None)

isovar_results = run_isovar(
    variants="cancer-mutations.vcf",
    alignment_file="tumor-rna.bam",
    read_collector=read_collector)
```

Unmapped reads, reads that failed vendor QC and records without a sequence are
always skipped. MAPQ is compared as a number, so STAR's 255 for a unique
alignment passes the default minimum of 1.
[Read eligibility across platforms](read-processing-audit.md) lists every rule,
the matching command-line options, and how quality tags from each platform are
treated.

## How reads are counted

- A read is one sequenced segment: one mate of a pair, or one long read. Its
  secondary and supplementary alignment records are not extra reads.
- A fragment is one sequenced template, so both mates of a pair count once.
  Fragment counts such as `num_alt_fragments` are usually the ones to filter on.
  Names are scoped by read group, and only complementary primary mates in the
  same read group are merged.
- A read whose alternative alignments support different alleles counts as
  other, not ref or alt, and cannot extend an assembly.
- None of these is a molecule count: Isovar does no UMI deduplication. The
  `*_read_names` properties are plain names for display.

## The RNA support record

Every support in the JSON exports (`protein-hypotheses`,
`allele-interpretations`, `sv-rna`, the SV ORF export and `fusion`) is one RNA
support record:

| Field | Meaning |
|---|---|
| `reads` | Sequenced segments: a mate, a single-end read or a long read |
| `fragments` | Templates: both mates of a pair count once, so this equals `reads` for single-end or long reads |
| `umis` | Distinct cell barcode and UMI pairs, within one library |
| `cells` | Distinct cell barcodes, within one library; a barcode without a UMI still counts |
| `umis_complete`, `cells_complete` | Whether there are reads and every one contributed, with a known library, so that the count is exact |
| `unlabeled_reads`, `unknown_library_reads` | Reads without a usable cell/UMI label, or whose library is unknown |
| `label_statuses` | Reads by label status |

The UMI and cell fields are null unless cell/UMI labels were assessed
(`--cell-umi-labels`; see [cell/UMI evidence](cell-umi-evidence.md)).

Tables use the same names, with a prefix where needed:

- `num_alt_reads`, `num_alt_umis` and `alt_umis_complete` in `allele-counts`;
- `protein_reads` and `translation_cells` in the `protein-hypotheses` TSV;
- bare `reads`, `umis` and so on in the SV ORF TSV.

`isovar plot` figure data count the read objects drawn, not this record.

## Overlapping mates, insertions and soft clips

Where overlapping mates disagree at a base, the higher-quality base wins. This
can change allele support as well as the assembled sequence. Mates stay
separate when their alignments conflict, or when a disagreement has equal or
missing qualities.

A read that ends at an insertion supports the reference allele only if both
flanking reference bases are aligned.

`use_soft_clipped_bases` keeps unaligned read ends; it does not realign them.

## Optional read annotations

These are off by default and leave raw counts unchanged:

- [Adapter and poly-A trimming](read-end-inference.md): the `ReadCollector`
  settings `infer_read_ends`, `read_end_profile`, `trim_adapters` and
  `trim_poly_a`. Original BAM records and aligned or inserted bases are never
  modified.
- [Base-quality assessments](base-quality.md) report focal-allele support that
  meets an explicit threshold alongside raw counts. Missing qualities remain
  unknown, and reconstruction uses the same pipeline across platforms.
- [Partial-read support scores](partial-read-support.md) divide each
  fragment's support among compatible RNA and protein hypotheses, and report
  per-cell summaries alongside raw counts.
- [Cell/UMI evidence](cell-umi-evidence.md) counts the cell barcodes and UMIs
  behind each allele and protein.
