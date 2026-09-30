# Cell/UMI evidence

Single-cell and UMI-tagged libraries label each read with a cell barcode (`CB`)
and a unique molecular identifier (UMI, such as `UB`). Several reads with the
same pair of labels may be copies of one RNA molecule. Isovar reports how many
distinct labels, and cells, support each piece of evidence, alongside the read
and fragment counts, under one policy, `isovar.cell_umi_labels.v1`:

- for small variants, each allele and protein hypothesis, and each allele
  interpretation (opt-in, below);
- for `reconstruct_sv_rna`, each junction, each path's voting reads and each
  ORF.

A label count is not a molecule count. Isovar uses the labels as the input
reports them. Label reporting does not cluster UMIs or correct barcodes. A missing or
unresolved label never changes ordinary sample reconstruction. The separate
opt-in [pool reconstruction](#pooled-reconstruction-with-per-cell-attribution)
analysis uses these labels to define reconstruction boundaries (#226).

## Small variants

```sh
isovar allele-counts --vcf variants.vcf --bam sc-rna.bam --cell-umi-labels \
    --sample-id tumor-1 --output counts.csv
isovar protein-hypotheses --vcf variants.vcf --bam sc-rna.bam --cell-umi-labels \
    --sample-id tumor-1 --output hypotheses.json
```

`allele-counts` adds these columns for each of `ref`, `alt` and `other`, beside
`num_*_reads` and `num_*_fragments`:

- `num_*_umis` and `num_*_cells`: distinct cell barcode and UMI pairs, and
  distinct cells;
- `*_umis_complete` and `*_cells_complete`: whether those counts are exact;
- `num_*_unlabeled_reads` and `num_*_unknown_library_reads`: reads without a
  usable label, and reads whose library is unknown;
- `num_cells_with_ref_and_alt`: cells with reads of both alleles.

`protein-hypotheses` fills in the cell/UMI fields of every support, alleles,
proteins and translations alike, and adds `cells_with_ref_and_alt` to each
event's `allele_support`. `allele-interpretations --cell-umi-labels` does the
same for each interpretation.

From Python, use `cell_umi_allele_evidence(results, alignment_file,
sample_id=..., source=...)`, or pass `cell_umi_alignment_file=` to
`export_protein_hypotheses`. Pass the `ReadCollector` the results were collected
with. A read whose record that collector would reject is reported as
`metadata_unavailable` and counted as unlabelled, with a warning.

A few cautions:

- Labels come from the eligible records overlapping the variant, so a mate
  outside the locus is not consulted for conflicts.
- Supports are counted separately, with no cell IDs, so do not add `cells`
  across alleles or proteins that may share cells.
- Counts from a selected or downsampled read set are not cell prevalence.
- If no read carries a usable `CB`, all cell counts are zero and a warning is
  logged; the data are probably not single-cell.

## Reading the counts

Every support, whether a small-variant allele or protein, an SV junction's
`direct_support` or an ORF's `full_interval_support`, is the
[RNA support record](../README.md#collecting-rna-reads). Its cell/UMI fields:

| Field | Meaning |
| --- | --- |
| `umis` | Distinct cell barcode and UMI pairs among the reads, within their declared library. Without `LB`, one pair in two read groups counts twice |
| `cells` | Distinct cell barcodes, by the same scoping; a barcode without a UMI still counts |
| `umis_complete` | Whether there are reads and every one has a usable pair and a known library, so `umis` is exact |
| `cells_complete` | Whether there are reads and every one has a trusted barcode and a known library, so `cells` is exact |
| `unlabeled_reads` | Reads without a usable complete label |
| `unknown_library_reads` | Reads without unambiguous library metadata |
| `label_statuses` | Reads by resolution status |

Use `umis` and `cells` as exact counts when their `*_complete` flag is true.
Otherwise the unlabelled and unknown-library counts show why they are not; with
no reads, both flags are false and both counts zero. A
label count is not a molecule count: UMIs are not clustered or corrected.

## Which labels count as the same

A label belongs to the input's `source` and `sample_id`. Within that input, the
read group header's `SM` and `LB` fields give the sample and library:

- Read groups with the same `SM` and `LB` share labels, so matching CB/UMI pairs
  in them count once.
- Different libraries or samples never share labels.
- Without `LB`, the library is unknown. Labels then stay local to their read
  group, or to the whole input when there is no `RG`. They count toward
  `umis` but leave `umis_complete` false.
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

In SV results, `cell_umi_evidence.reads` lists every read consulted, with its scope,
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

## Pooled reconstruction with per-cell attribution

`protein-hypotheses --cell-reconstruction sample` adds a separate
`cell_reconstruction` object to the JSON. The ordinary sample result and TSV
keep their existing meaning. This opt-in analysis reconstructs local protein
candidates within explicit pools and compares the original observations with
the discovered candidate catalog:

```sh
isovar protein-hypotheses --vcf variants.vcf --bam sc-rna.bam \
    --sample-id tumor-T1 --source library-set-T1 \
    --cell-reconstruction sample --cell-max-edits 1 \
    --output hypotheses.json
```

Python callers can use `reconstruct_cell_groups(results, alignment_file,
sample_id=..., source=..., mode="sample")`. Pass the original `ReadCollector`
and transcript whitelist if customized. The input results must record their
`protein_sequence_settings`. Discovery reuses those settings but keeps all
protein hypotheses, independently of the ordinary output's top-protein cap.
The configured per-base coverage floor applies inside each pool; it is never
lowered to obtain a per-cell result.

The pooling modes are:

| Mode | Discovery scope |
| --- | --- |
| `sample` | Explicit sample/source plus header `SM`; different declared samples stay separate |
| `library` | Existing sample/library namespace, with RG fallback when the library is unknown |
| `cell` | One trusted, scoped barcode; unresolved cells are excluded from discovery |
| `group` | An external cell-group name within a declared sample |

For group mode, supply `--cell-groups groups.json`, a JSON object mapping the
exported `cell_id` values to group names. Obtain those IDs from a first sample
export with the same `--sample-id` and `--source`. Identical barcode strings
from different libraries have different IDs. Explicit group membership may
pool their reads for discovery, but it never merges their cell identities.
Cells absent from the mapping are excluded from discovery and still receive
attribution against candidates from their declared sample. Supply distinct
sample/source identities for timepoints and separate biological inputs; missing
metadata cannot recover distinctions that were already lost upstream.

Each event reports reconstruction pools (including no-protein reasons), the
union of nucleotide candidates and their protein hypothesis IDs, and a row for
each observed trusted cell. Reads without trusted cell identities are included
in `unresolved_cell_support`; reads without the requested pool are included in
`excluded_from_pool`. These sets may overlap. No read is silently dropped.
The export includes hashed evidence sets for tracing and deduplicating support.

Attribution is separate from discovery. Four cells with one informative read
each can reconstruct a candidate at a two-read coverage floor and each retain
one supporting observation. A noisy read excluded from assembly can still be
compatible during attribution. Conversely, reads covering only shared sequence
remain ambiguous between candidates. Cells can support multiple sequences.

`score_policy=anchored_edit_compatibility.v1` compares the observed prefix,
allele and suffix with each candidate, anchored at the nominated variant.
Only overlapping flanks are compared. Each part uses unit-cost edit distance;
splice compatibility is checked against candidate transcripts. The default
maximum total cost is 1 (`--cell-max-edits`); the default minimum number of
compared read bases is 10 (`--cell-min-overlap`). Non-ACGT observations in the
compared interval remain unassessed. The `comparisons` report costs and the
candidate interval covered; `spans_candidate` refers to that local nucleotide
candidate, potentially supported by merged mates, not necessarily a full ORF.

Every candidate also has a paired reference-allele competitor with the same
flanks. An observation whose reference cost is equal or better cannot provide
mutant-specific support. All candidates within the tolerance remain possible;
the implementation does not force a winner by majority abundance. Alternative
placements of one observation must agree for a unique assignment.

`unique_within_catalog` means one compatible candidate after these checks.
`ambiguous` retains the candidate set, including disagreements between
placements. Other statuses report no candidate, inadequate coverage/path,
reference-or-other allele evidence, or excessive edit distance. Candidate
`unique_cell_ids` and `ambiguous_cell_ids` list cells with at least one such
observation; a cell may appear in both lists or support multiple candidates.
Do not add these lists or interpret them as mutually exclusive cell states.
`num_unique_cells` and `num_ambiguous_cells` are the corresponding observed
counts. `protein_attribution` aggregates the same evidence by protein hypothesis:
ambiguity between synonymous nucleotide candidates can still identify one
protein. Observation `hypothesis_ids` and `protein_status` record this decision.
Multiple protein hypotheses for one nucleotide candidate remain unresolved.

These are compatibility costs, not posterior probabilities or validated
expression calls. The first implementation does not model base qualities,
ambient RNA, barcode errors or systematic platform errors. Reference competition
covers the focal allele in each discovered context, not an exhaustive catalog
of wild-type isoforms. Discovery still uses the existing exact-overlap assembler;
edit-tolerant attribution does not perform noisy-read consensus correction.
Cells without an informative read are unknown, and cells absent from the locus
are not enumerated. Observed cell counts are not an estimate of the total number
of expressing cells. Local RNA support does not establish protein translation.

The discovery/quantification separation follows
[IsoQuant](https://ablab.github.io/IsoQuant/cmd.html) and
[FLAMES](https://doi.org/10.1186/s13059-021-02525-6).
[Longcell](https://doi.org/10.1038/s41467-025-60902-2) describes the UMI and mapping
errors that require additional modeling for quantitative single-cell inference.
