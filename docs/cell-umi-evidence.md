# Cell/UMI evidence in SV reconstruction

`reconstruct_sv_rna` reports cell/UMI **labels**, separate from sequenced
segments, templates, ONT signal ancestry and independent molecules, under the
policy `isovar.cell_umi_labels.v1`. Junction and ORF support use the same
accounting. It does not perform cell-stratified reconstruction
([#226](https://github.com/openvax/isovar/issues/226)), UMI clustering or
barcode correction.

## Scope

Every label belongs to the explicit input `source` and `sample_id`. Within that
namespace, header `SM` and `LB` identify the declared sample and library. Two
read groups with the same `SM`/`LB` share a label if their complete CB/UMI pairs
match. Different libraries or header samples never merge. When `SM` is absent,
the explicit input sample still supplies the sample namespace; missing `SM`
is retained and is not equated with a different group's populated `SM`.

Missing `LB` means unknown library scope. Labels then remain read-group-local,
or input-local when `RG` is absent. They contribute to `observed_labels`, but
not a complete library-scoped count. Duplicate header RG IDs are ambiguous and
their labels remain unresolved. No donor, timepoint, library or barcode suffix
is inferred from a filename. Header metadata are declarations, not independently
validated biological sample assignments.

Sources remain separate even when their label strings agree. Different input
paths do not prove independent evidence: do not add counts from reprocessed
products or otherwise overlapping read sets. Pooled donor identity, molecule
collisions and cross-source reconciliation remain unassessed.

## Tags and conflicts

- Use `CB` and `UB` verbatim as reported cell/UMI identifiers. Do not strip
  suffixes, change case, correct strings or substitute raw `CR`/`XC`, `UR`/`OX`
  or possibly raw `RX`. SAM's CB is *optionally* corrected; the presence of a
  tag alone does not prove whitelist correction or cell calling. `UB` uses
  the documented corrected-UMI convention without claiming its producer has
  been independently validated. Gene-dependent UMI correction is not rerun.
- `XM` contributes only when a record's PG pointer, its RG's PG pointer, or an
  unambiguous whole-header program history identifies an Iso-Seq `correct`
  step. Follow `PP` links; reject missing, cyclic, contradictory and duplicate
  program references. Without a pointer, every header lineage must agree.
  A later Iso-Seq `tag` step makes XM raw; `dedup`/`groupdedup` alone does not
  establish correction. Platform labels, filenames and an unrelated Iso-Seq
  header entry are insufficient. Bismark's XM is methylation, not a UMI.
- All retained placements of a segment participate in conflict checking.
  Non-string/empty selected tags and disagreeing CB or accepted UMI values
  leave the segment unresolved. Corrected UB and attributable corrected XM
  must agree when both occur. Raw and corrected values may differ normally.
- At least one record must carry a complete pair. CB on one record and UMI
  on another do not become a complete label. Missing tags on additional
  records are counted explicitly. Conflicting labels on visible mates leave
  the template's segments unresolved, even if only one mate is a witness.

The original records and native evidence tags remain unchanged. A missing or
unresolved label never removes an RNA alignment or changes reconstruction.

## Counts

Both `junctions[].direct_cell_umi_support` and
`exploratory_orfs.candidates[].full_interval_support.cell_umi_support` contain:

| Field | Meaning |
| --- | --- |
| `unit` | `cell_umi_label` |
| `segment_ids` | Exact supporting `(RG, QNAME, mate bits)` identities |
| `observed_labels` | Distinct resolved label pairs in their declared scopes; a subset count |
| `unresolved_segments` | Witnesses without a usable complete label |
| `unknown_library_segments` | Witnesses without unambiguous declared library metadata |
| `all_segments_labeled` | Whether every witness has a usable label, independent of library completeness |
| `complete_label_count` | Count only if there are witnesses and every label and library scope is resolved; otherwise null |
| `independent_molecules` | Always null; not established by this policy |
| `status_counts` | Segment counts for each resolution status |

`cell_umi_evidence.segments` records consulted segment identities, scope,
selected labels, UMI tags, XM semantics and reasons. Its label key is the tuple
`(source, sample_id, header_sample, library, read_group, basis,
cell_barcode, umi)`, where `read_group` is set only when the sample/library
scope is unknown and `basis` is `sample_library`, `read_group` or `input`. Unknown metadata are null. Visible mates consulted for
conflicts can appear here without appearing in a support denominator.

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
