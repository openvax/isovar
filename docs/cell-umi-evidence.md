# Cell/UMI evidence

Single-cell and UMI-tagged libraries label each read with a cell barcode (`CB`)
and a unique molecular identifier (UMI, such as `UB`). Several reads with the
same pair of labels may be copies of one RNA molecule. Isovar reports how many
distinct labels, and cells, support each piece of evidence, alongside the read
and fragment counts, under one policy, `isovar.cell_umi_labels.v1`:

- for small variants, each allele and protein hypothesis (opt-in, below);
- for `reconstruct_sv_rna`, each junction and each ORF.

A label count is not a molecule count. Isovar uses the labels as the input
reports them. It does not cluster UMIs, correct barcodes or reconstruct each
cell separately ([#226](https://github.com/openvax/isovar/issues/226)). A
missing or unresolved label never removes a read or changes the reconstruction.

## Small variants

```sh
isovar allele-counts --vcf variants.vcf --bam sc-rna.bam --cell-umi-labels \
    --sample-id tumor-1 --output counts.csv
isovar protein-hypotheses --vcf variants.vcf --bam sc-rna.bam --cell-umi-labels \
    --sample-id tumor-1 --output hypotheses.json
```

`allele-counts` adds these columns for each of `ref`, `alt` and `other`:

- `num_*_cells` and `num_*_cell_umi_labels`: distinct cells and cell/UMI labels
  (lower bounds);
- `num_*_unlabeled_segments`: reads without a usable label;
- `num_cells_with_ref_and_alt`: cells with labelled reads of both alleles.

`protein-hypotheses` adds an event-level `cell_umi_evidence`, and a
`cell_umi_support` in every protein's `rna_support`.

From Python, use `cell_umi_allele_evidence(results, alignment_file,
sample_id=..., source=...)`, or pass `cell_umi_alignment_file=` to
`export_protein_hypotheses`. Each summary below also gives `observed_cells`,
and `complete_cell_count` under the same rule as `complete_label_count`.
Supports carry an `evidence_set_id` rather than read names.

Labels come from the eligible records overlapping the variant. A mate outside
the locus is not consulted for conflicts. Counts from a selected or downsampled
read set are not cell prevalence.

## Reading the counts

Each summary (a small-variant allele or protein, an SV junction's
`direct_cell_umi_support`, or an ORF's `full_interval_support.cell_umi_support`)
contains:

| Field | Meaning |
| --- | --- |
| `complete_label_count` | Distinct labels among the supporting reads, if every read has a usable label and a known library; otherwise null |
| `observed_labels` | Distinct usable labels found, even when some reads lack one; a lower bound |
| `all_segments_labeled` | Whether every supporting read has a usable label, whether or not its library is known |
| `unresolved_segments` | Supporting reads without a usable complete label |
| `unknown_library_segments` | Supporting reads without unambiguous library metadata |
| `status_counts` | Supporting reads by resolution status |
| `segment_ids` | The supporting reads, as `(RG, QNAME, mate bits)`. SV results only; exports and small-variant summaries use hashed IDs instead |
| `unit` | Always `cell_umi_label` |
| `independent_molecules` | Always null: labels do not establish independent molecules |

Use `complete_label_count` when it is not null. When it is null, look at
`observed_labels` and the unresolved and unknown-library counts to see why.

## Which labels count as the same

A label belongs to the input's `source` and `sample_id`. Within that input, the
read group header's `SM` and `LB` fields give the sample and library:

- Read groups with the same `SM` and `LB` share labels, so matching CB/UMI pairs
  in them count once.
- Different libraries or samples never share labels.
- Without `LB`, the library is unknown. Labels then stay local to their read
  group, or to the whole input when there is no `RG`. They count toward
  `observed_labels` but never toward `complete_label_count`.
- Without `SM`, the input's `sample_id` names the sample. A missing `SM` is not
  treated as equal to another group's `SM`.
- Duplicate read group IDs in the header are ambiguous, and their labels stay
  unresolved.
- Nothing is inferred from file names: no donor, timepoint, library or barcode
  suffix.

Header fields are declarations, not independently validated sample assignments.
Different inputs are never combined, even when their label strings agree.
Different input files are not independent evidence either: do not add counts
from reprocessed products or overlapping read sets. Pooled donors, UMI
collisions and reconciliation across inputs are not assessed.

## Which tags are used

- `CB` and `UB` are used verbatim: suffixes are not stripped, case is not
  changed and strings are not corrected. The raw tags (`CR`/`XC`, `UR`/`OX`, and
  possibly `RX`) are not substituted. In the SAM specification CB is only
  *optionally* corrected, so the tag alone does not prove whitelist correction
  or cell calling. `UB` is taken as the documented corrected UMI, without
  independent validation of its producer. Gene-dependent UMI correction is not
  rerun.
- `XM` counts only when the BAM header shows it was written by an Iso-Seq
  `correct` step. That step must be named by the record's `PG` pointer, its read
  group's `PG` pointer, or an unambiguous program history for the whole header.
  - `PP` links are followed. Missing, cyclic, contradictory and duplicate
    program references are rejected.
  - Without a pointer, every program lineage in the header must agree.
  - A later Iso-Seq `tag` step makes `XM` raw again. `dedup` or `groupdedup`
    alone does not establish correction.
  - Platform labels, file names and an unrelated Iso-Seq header entry are not
    enough. Bismark's `XM` is a methylation call, not a UMI.

## Conflicts

A read counts only if at least one of its records carries both a cell barcode
and a UMI. A CB on one record and a UMI on another do not make a complete label.
Records missing tags are counted explicitly. A read is left unresolved when:

- any retained placement has an empty or non-string tag;
- its placements disagree on CB or on the accepted UMI;
- corrected `UB` and attributable corrected `XM` disagree (raw and corrected
  values may differ normally);
- its visible mates carry conflicting labels, even if only one mate supports
  the junction.

## The evidence record

`cell_umi_evidence.segments` lists every read consulted, with its scope,
selected labels, UMI tags, `XM` interpretation and the reasons for its status.
Visible mates checked for conflicts can appear here without being counted as
support. Each label's key is the tuple `(source, sample_id, header_sample,
library, read_group, basis, cell_barcode, umi)`:

- `basis` is `sample_library`, `read_group` or `input`, the scope the label was
  resolved in;
- `read_group` is set only when the sample or library is unknown;
- unknown metadata are null.

The original records and their tags are not modified.

## Primary sources

- [SAM header RG/SM/LB definitions](https://samtools.github.io/hts-specs/SAMv1.pdf)
  and [barcode tags](https://samtools.github.io/hts-specs/SAMtags.pdf).
- [Cell Ranger BAM tags](https://www.10xgenomics.com/support/software/cell-ranger/latest/tutorials/outputs/cr-outputs-bam):
  CB suffixes and the cell/gene context of UB correction.
- [Iso-Seq tags](https://isoseq.how/isoseq-tags.html) and
  [deduplication](https://isoseq.how/umi/dedup-faq.html): correction stages and
  the sequence/alignment evidence used beyond matching barcode strings.
- [Bismark alignment output](https://github.com/FelixKrueger/Bismark/blob/master/docs/src/content/docs/usage/alignment.md):
  XM methylation-call semantics.
