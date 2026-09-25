# Changelog

Behavior changes that can alter results or break callers, by release. Patch
releases that only fix bugs or add fixtures are omitted; see the
[commit history](https://github.com/openvax/isovar/commits/master) for those.

## 1.37.0

Every output reports supporting RNA with one record, and "segments" is renamed
"reads" throughout.

- **The record.** It has `reads`, `fragments`, `umis`, `cells`,
  `umis_complete`, `cells_complete`, `unlabeled_reads`, `unknown_library_reads`
  and `label_statuses`. The UMI and cell fields are null unless cell/UMI labels
  were assessed. It is defined once, in `isovar.rna_evidence.rna_support`, and
  described in the README.
- **Where it replaces the earlier forms:**
  - SV junctions: `direct_support` replaces `direct_segments`,
    `direct_fragments` and `direct_cell_umi_support`.
  - SV paths: `sequence_evidence.voting_support` replaces the voting counts.
  - SV ORFs: `full_interval_support`, and splice-inclusion `support`, carry
    the record's fields.
  - Fusions: the junction has `direct_support`, and `evidence.support` covers
    all reads.
  - Protein hypothesis exports: `rna_support` and every `allele_support`
    entry are records.
  - Allele interpretation supports are records too.
- **Removed:** `observed_labels`, `complete_label_count`, `independent_molecules`,
  `cell_umi_support`, the export's `cell_umi_alleles`, and raw `segment_ids`
  in summaries.
- **Renamed:** lineage summaries use `reads`, `signal_groups` and
  `unresolved_reads`. Evidence sets hold `read_ids`. `*_segments` counts become
  `*_reads` (`missing_quality_reads`, `read_path_notes`, and so on).
- **Schema versions:** `sv_rna_candidates.v5`, `sv_rna_orfs.v4` (with new TSV
  count columns), `fusion_rna.v3`, `protein_hypotheses.v2` and
  `allele_interpretations.v2`.
- **`isovar allele-counts --cell-umi-labels`** columns are
  `num_*_umis`, `num_*_cells`, `num_*_unlabeled_reads` and
  `num_*_unknown_library_reads`.
- **`isovar allele-interpretations --cell-umi-labels`** counts UMIs and cells.
- **Varcode:** requires `>=10,<11`. The suite passes on 10.3.0.

## 1.36.0

- `IsovarReadPhasing.in_cis(v1, v2)` reports cis or trans from fragments that
  cover both loci, and `None` otherwise. Trans needs fragments carrying one alt
  allele with the other variant's reference allele. A matched germline variant
  is cis when its edit is in the variant's top assembled protein. Varcode's
  `MolecularPhaseResolver` calls this method instead of reading trans from the
  partner lists, which only compare variants in the run
  ([varcode#517](https://github.com/openvax/varcode/issues/517)).

## 1.35.0

- Small variants count the cells and cell/UMI labels behind their allele reads
  and protein hypotheses. This is the first stage of
  [#226](https://github.com/openvax/isovar/issues/226), and it is opt-in:
  sample-level reconstruction is unchanged.
  - Python: `cell_umi_allele_evidence` and `isovar.cell_evidence.CellUmiAlleles`.
  - `isovar allele-counts --cell-umi-labels --sample-id` adds these columns for
    each allele: `num_*_cells`, `num_*_cell_umi_labels`,
    `num_*_unlabeled_segments` and `num_*_unknown_library_segments`. It also
    adds `num_cells_with_ref_and_alt`.
  - `isovar protein-hypotheses --cell-umi-labels` (`cell_umi_alignment_file=`
    in Python) adds `cell_umi_alleles` to each event and `cell_umi_support` to
    each protein's `rna_support`.

  Labels use the SV policy `isovar.cell_umi_labels.v1`. Cells are counted
  from a trusted barcode, with or without a UMI. Reads whose record the
  collector would reject are unresolved (`metadata_unavailable`) and logged.
- The cell/UMI docs no longer call `observed_labels` a lower bound when
  libraries are unknown: without `LB`, one label in two read groups counts
  twice.

## 1.34.0

- `reconcile_allele_interpretations` and `isovar allele-interpretations` compare
  competing interpretations of one locus with the RNA reads, as
  `isovar.allele_interpretations.v1`. This is the first stage of
  [#306](https://github.com/openvax/isovar/issues/306). An interpretation can
  be an SNV, a multi-base change, several adjacent variants or a deletion.
  Each is `supported`, `contradicted_by_informative_evidence`,
  `insufficient_RNA` or `unassessed`. Interpretations with the same RNA
  sequence are reported as indistinguishable. Every observed allele is listed
  with its difference from the reference, and support uses hashed evidence
  sets. On osteosarc RNA, NTF3's reads carry Mutect2's AG>GT (the same as
  Strelka's two SNVs), not the catalogued A>G. GLIS3's reads carry the
  catalogued deletion together with a somatic A>T.

## 1.33.0

Assembled cDNA edits are attributed to supplied variants
([#297](https://github.com/openvax/isovar/issues/297)):

- `run_isovar(germline_variants=...)` and `--germline-vcf` take variants from a
  matched normal. An edit in an assembled cDNA that is one of them is now in
  `ProteinSequence.known_germline_transcript_edits`.
- An edit that is another input variant is co-somatic and is now in
  `known_somatic_transcript_edits`, where before it was unexplained. Only
  edits that match no supplied variant stay in `unexplained_transcript_edits`.
- Matching happens on the transcript. It accepts indels shifted within a repeat,
  and splits runs of adjacent changes into the supplied variants they contain.
  A variant may not be both somatic and germline; an edit that is equivalent
  to one of each stays unexplained.
- `ProteinSequence` gains `known_variants` and `with_known_variants`.
  `IsovarResult` gains the co-somatic and germline variant sets of the top
  protein, and `num_*` counts of them and of unexplained edits, which are new
  `isovar run` columns and filters.
- `PhaseGroup.germline_variants` lists the germline variants in the group's
  assemblies. Protein hypothesis exports label each edit `nominated_variant`,
  `co_somatic_variant`, `germline_variant` or `unexplained`, with its
  `source_event_id`.
- `IsovarMutantTranscript` evidence counts reads and fragments by read-group
  identity, like other counts, instead of by name.
- A `min_`/`max_` filter on a property that is None, such as a top-protein
  measure for a variant with no protein, now fails instead of raising
  `TypeError`.

## 1.32.0

- `export_protein_hypotheses` and `isovar protein-hypotheses` export every
  protein Isovar reconstructs for each variant, as `isovar.protein_hypotheses.v1`
  JSON and TSV ([#324](https://github.com/openvax/isovar/issues/324)). Each
  protein keeps all its translations, including synonymous cDNAs, and marks the
  shorter windows another protein contains. RNA support is given as counts and
  an evidence set of hashed read identities. The command keeps every protein by
  default.
- `union_rna_support` combines read sets from several hypotheses or exports,
  counting each read once, and refuses sets from different evidence scopes.
  `isovar sv-rna` ORF exports use the same read IDs, and their `rna_support`
  gains `evidence_scope`.
- `ProteinSequenceCreator.settings()` returns the creator's effective
  parameters. `run_isovar` records them on each result as
  `IsovarResult.protein_sequence_settings`.

## 1.31.2

- CI publishes coverage with the `coveralls` package it already installs from
  PyPI, instead of an action that downloads a separate reporter binary; that
  download intermittently failed checksum verification
  ([#358](https://github.com/openvax/isovar/issues/358)). A failed upload still
  fails the job.

## 1.31.1

- The README and guides are rewritten for readers. Each guide opens with what it
  is for and a quick start, and dense paragraphs became tables and lists.
  `docs/read-processing-audit.md` is now a reference on which reads are used, and
  the README groups the guides by task, and explains how to choose the reference
  annotation from Python ([#123](https://github.com/openvax/isovar/issues/123)).
  No behavior changes.
- The docs test also checks that `#section` links point to an existing heading.

## 1.31.0

Supplied fusions and reconstructed SV paths share one RNA path format
([#364](https://github.com/openvax/isovar/issues/364)). Neither output had
production consumers, so the old forms are not kept:

- `reconstruct_fusion` returns schema `isovar.fusion_rna.v2`. Its supplied
  transcript is `paths[0]`, with the same `sequence`, `junctions`, `frame_status`
  and `translations` fields as `reconstruct_sv_rna`. Renamed:
  - `schema_version` → `schema`;
  - `cdna_sequence` → `sequence`;
  - `junction_interval` → the junction's `query_interval`;
  - `junction_peptides` → `candidate_peptides`;
  - `donor_transcript_ids` → `transcript_ids`;
  - `directly_spanning_fragments` → `direct_fragments`;
  - status `unresolved_frame` → `unresolved`.

  Translations gain `translation_end` and `translation_observed`.
- `reconstruct_sv_rna` returns schema `isovar.sv_rna_candidates.v4`:
  - junction and inclusion-edge `query_interval`s are half-open intervals of
    the unplaced bases between partners, instead of the two flanking offsets;
  - translations gain a top-level `complete_5prime`;
  - `reference_comparisons[].start_kind` uses the `start_evidence` region
    names, so `annotated_start` becomes `annotated_CDS_start`;
  - `RnaObservation.reverse` becomes `reverse_complement`;
  - junction `direct_molecules` and ORF `molecule_labels` are removed. They
    repeated `complete_label_count`, a count of cell/UMI labels rather than
    molecules; read it from `direct_cell_umi_support` or `cell_umi_support`.
- ORF exports are `isovar.sv_rna_orfs.v3` with one `interval_convention`, and
  prediction comparisons are `isovar.sv_rna_prediction_comparison.v2`.
  `export_sv_rna_orfs` accepts only v4 reconstructions, and `write_sv_rna_orfs`
  only v3 exports. `normalize_sv_rna_orf_export` and its v1 flag migration are removed.
- `fusion_from_dict` and `sv_rna_input_from_dict` reject unknown top-level keys.

## 1.30.0

- `isovar sv-rna`'s splice-inclusion gate treats MAPQ 255 as a unique alignment,
  as STAR writes it, so STAR-aligned intronic starts can reach priority 3
  ([#361](https://github.com/openvax/isovar/issues/361)). `inclusion_mapq_255_is_unique`
  now defaults to True. `--inclusion-mapq-255-unavailable` restores the SAM
  specification's meaning and replaces `--inclusion-mapq-255-is-unique`.

## 1.29.2

- CI rejects a pull request whose version is not newer than its base branch's, or
  is already on PyPI or tagged ([#345](https://github.com/openvax/isovar/issues/345)).

## 1.29.1

- The complete osteosarc figure gallery's README names its own generator
  (`examples.osteosarc_context_figures`) and links its reports; the base-only
  command is labelled as such ([#301](https://github.com/openvax/isovar/issues/301)).

## 1.29.0

Consolidated SV helpers ([#364](https://github.com/openvax/isovar/issues/364));
reconstruction output is unchanged.

- `sv_rna_orfs.exploratory_orfs` takes the supplied `FusionReference` models
  instead of the private `sv_rna._Model` wrappers.
- `compare_sv_rna_predictions` is exported from `isovar`.
- `isovar.sv_rna_relations` defines the junction linkage statuses and the
  event-linked and event-crossing subsets once. They were previously spelled out
  separately in `sv_rna`, `sv_rna_orfs` and `sv_rna_comparison`.
- `chimeric_alignment.sa_entries` is the one `SA` tag parser. SV record retrieval
  and chimeric-path validation both use it. Path validation no longer
  suppresses `AttributeError`/`TypeError`, which hid programming errors.
- `segment_identity` moved to `read_identity`; `sv_rna` still imports it.
- `read_metadata.ProgramHistory` walks `@PG` chains and resolves record/read-group
  `PG` pointers for both the cell/UMI and ONT lineage evidence. Each keeps its own
  producer policy.
- Each SV record's read-end view is built once rather than twice, which halves
  adapter matching when a read-end profile is set.

## 1.28.1

- Package metadata declares the license as the SPDX expression `Apache-2.0`
  (PEP 639) instead of the deprecated classifier; building requires setuptools 77
  or later ([#205](https://github.com/openvax/isovar/issues/205)).

## 1.28.0

CSV output from the command-line tools
([#258](https://github.com/openvax/isovar/issues/258)):

- Every table starts with `variant`, `chr`, `pos`, `ref`, `alt`. `chr` keeps the
  input's contig name (`chr9`, not `9`), and `isovar run` gains the four columns
  after `variant`. Existing `isovar run` columns are unchanged.
- Collection fields were written as counts under their original names. They are
  now `num_translations` (protein sequences), `num_reads` and `num_fragments`
  (variant sequences), and `compatible_transcript_ids` lists IDs.
- `isovar translations` renames `variant_orf` to `in_frame_cdna_sequence` and
  `reference_context` to `reference_transcript_names`, matching their contents.
- `isovar protein-sequences` adds `num_supporting_reads`, `num_supporting_fragments`
  and `transcript_ids`.
- Floats are written with six significant digits.
- An empty result keeps its header row; `isovar_results_to_dataframe([])` returns
  the `IsovarResult.record_columns()` columns.

## 1.27.1

- Command-line input problems are one-line usage errors (exit 2) instead of
  tracebacks ([#256](https://github.com/openvax/isovar/issues/256)): missing
  variant/BAM/JSON files, missing `--genome`, SAM or unindexed alignment files,
  a missing output directory (checked before any work), unknown
  `--output-columns`, invalid `ReadCollector` settings and out-of-range numeric
  options (negative lengths or counts, fractions outside [0, 1], `--dpi` below 72).
  `isovar.cli.validation.CommandInputError` is a `ValueError`.
- `fusion_from_dict` and `sv_rna_input_from_dict` raise `ValueError` for a
  missing key or unexpected field instead of `KeyError`/`TypeError`, so
  `isovar fusion` and `isovar sv-rna` no longer report programming errors as
  usage errors.

## 1.27.0

- `isovar sv-rna` output schema `isovar.sv_rna_candidates.v3`
  ([#363](https://github.com/openvax/isovar/issues/363)). Changes since v2:
  - each path has `end_reasons` (`observations_end`, `repeated_genomic_position`
    or `path_limit` for each end), so a truncated path no longer looks complete;
  - `parameters.read_collection` records every read-collection setting,
    including duplicate/secondary/QUAL policy, trimming and read-end profile
    SHA-256s; `use_soft_clipped_bases`, `custom_read_filter` and
    `min_mapping_quality` moved there from the top level of `parameters`;
  - `parameters.inclusion` and the 1.26.0 inclusion `thresholds` changes;
  - fields added after v2 was introduced: `start_evidence`, `junction_crossings`,
    `reconstruction_scopes` and `junction_links`, and the 1.25.1
    `departs_elsewhere_only` frame status;
  - `direct_molecules` and ORF `molecule_labels` are null unless every label and
    library scope is resolved (since 1.21.9);
  - records without a read name are reported as `unnamed_records_skipped`.
- `export_sv_rna_orfs` accepts v2 and v3 results. Occurrences carry `end_reasons`,
  and a truncated path end adds the `path_end_truncated` warning.

## 1.26.0

- `isovar sv-rna` splice-linked inclusion gates are configurable and recorded:
  `--inclusion-min-splice-anchor-bases` (8), `--inclusion-min-base-quality` (20),
  `--inclusion-min-mapping-quality` (20) and `--inclusion-mapq-255-is-unique`,
  also `reconstruct_sv_rna(inclusion_*=...)`. The MAPQ gate previously ignored
  every setting and was hard-coded at 20. Defaults are unchanged: MAPQ 255 is
  unavailable, and excluding it adds `mapq_255_excluded_from_inclusion` to
  `limitations`. For STAR input, the new flag lets unique alignments reach
  priority 3 ([#361](https://github.com/openvax/isovar/issues/361)).
- `annotate_orf_inclusion`'s `min_anchor_bases` is renamed
  `min_splice_anchor_bases`, as is the matching `thresholds` key, so it is not
  confused with the 18-base frame anchor. Witness rows gain `mapping_quality_255`,
  and their `minimum_mapping_quality` is the lowest MAPQ other than 255.

## 1.25.1

- `isovar sv-rna` ORF starts: only intronic start assessments carry
  `splice_inclusion` and `splice_inference="observed_path_assessed"`. UTR, CDS,
  noncoding and antisense starts were reported as assessed-but-unresolved.
- `frame_status` is `departs_elsewhere_only`, not `no_gene_anchor`, when a path has
  an exact coding anchor but every reading leaves the model before a reported
  junction ([#362](https://github.com/openvax/isovar/issues/362)).

## 1.25.0

- Importing Isovar no longer configures logging. Previously every module ran
  `logging.config.fileConfig`, which disabled loggers the host application had
  already created and sent Isovar, Varcode and PyEnsembl INFO logs to stdout.
  Commands now log to stderr, accept `--log-level`, and report per-candidate
  detail only at DEBUG. The unused `isovar/logging.conf` is removed.
- `ProteinSequenceCreator` instances with different settings no longer compare
  equal. `sort_protein_sequences` defaults to the pipeline's `balanced`
  preference rather than `support`.
- `--min-protein-sequence-support-fraction` belongs to the protein-sequence
  commands; `isovar translations` never ranks proteins and ignored it.
- `isovar allele-reads` writes `isovar-allele-reads-result.csv` by default
  instead of `output.csv`. `--help` shows every command's default output file.
- `isovar fusion` writes its result before plotting, so a plotting error keeps it.
- `reverse_complement_dna` complements N and IUPAC ambiguity codes in either case
  instead of raising `KeyError`.
- `sv-rna` no longer raises `KeyError('ID')` for `@RG`/`@PG` header lines without an ID.
- Removed: `filter_variant_sequences`, `filter_variant_sequences_by_length` and
  `filter_variant_sequences_by_read_support` (use `trim_variant_sequences`); the
  phasing helpers `create_read_names_to_variants_dict`,
  `create_variant_to_alt_read_names_dict`,
  `create_variant_to_protein_sequence_read_names_dict`, `compute_phasing_counts`,
  `threshold_phased_variant_counts` and `create_phase_groups`; the CLI helpers
  `read_evidence_dataframe_from_args` and `variants_reads_dataframe_from_args`;
  `read_end_inference.reverse_complement` (use `dna.reverse_complement_dna`).
- `PROTEIN_SELECTION.md` moved to `docs/protein-selection.md`.

## 1.21

- 1.21.7: Small-variant collection rejects vendor-QC-failed records (SAM flag
  0x200), as SV collection already did. `ReadCollector(read_filter=...)` applies
  a caller-supplied record predicate in both pipelines.
- 1.21.0: `isovar sv-rna` output schema `isovar.sv_rna_candidates.v2`: v1's
  `spanning_*` fields became `direct_*`, and whole-path containment counts became
  voting counts. Reconstruction scales to deep long-read loci; on Sid ONT T1
  TPST1–CRCP (8k records) 1.20.0 took 40 s and 1.9 GB, 1.21 takes 5 s and 0.3 GB.

## 1.20.0

- `isovar sv-rna` / `reconstruct_sv_rna` reconstructs RNA paths around one
  nominated SV ([guide](docs/sv-rna.md)).

## 1.19.0

- Opt-in adapter/poly-A inference and trimming (`--infer-read-ends`,
  `--read-end-profile`, `--trim-adapters`, `--trim-poly-a`). Defaults unchanged
  ([guide](docs/read-end-inference.md)).

## 1.18

- 1.18.1: A read ending at an insertion supports the reference allele only when
  both flanking reference bases are aligned.
- 1.18.0: Reads without QUAL are kept, with unknown base quality, instead of being
  discarded. `--require-base-qualities` /
  `ReadCollector(use_reads_without_base_qualities=False)` restores the old policy.

## 1.17

Read support counts sequenced segments, not alternative SAM alignments.

- 1.17.3: Supplementary pieces of one segment phase when reciprocal `SA` tags
  declare the same chimeric path ([#286](https://github.com/openvax/isovar/issues/286)).
- 1.17.2: Phasing requires a compatible pair of placements for a shared fragment
  ([#284](https://github.com/openvax/isovar/issues/284)).
- 1.17.1: Phasing uses read-group-scoped fragment IDs
  ([#282](https://github.com/openvax/isovar/issues/282)).
- 1.17.0: Only complementary primary mates in one read group are collapsed;
  conflicting allele calls from one segment become `other` reads
  ([#264](https://github.com/openvax/isovar/issues/264)). This can change counts,
  filtering and reconstructed context for multimapped reads.

## 1.16.0

- BAM-derived RNA prefixes carry the coding frame through their observed alignment,
  so an upstream indel is no longer hidden by equal-length trimming
  ([#265](https://github.com/openvax/isovar/issues/265),
  [details](docs/aligned-reading-frame.md)).

## 1.15.0

- `isovar fusion` / `reconstruct_fusion` validates supplied fusion transcripts
  ([guide](docs/fusion.md)).
- 1.14.1: Symbolic SVs and breakends are rejected by the small-variant pipeline.

## 1.11.0

- `ReadCollector()` and CLI collectors merge overlapping mates by default, matching
  `run_isovar()`; `--no-merge-overlapping-fragments` keeps them separate.
- The CLI applies the complete default filter set used by `run_isovar()`.
- `--reference-context-size` is accepted only by `isovar reference-contexts`.

## 1.10.0

- Reads keep their aligned exon blocks and splice junctions; assembly and
  translation are restricted to compatible transcripts.

## 1.9.0

- Overlap assembly is on by default in the Python API and CLI. Pass
  `variant_sequence_assembly=False` or `--disable-variant-sequence-assembly` to disable it.

## 1.8.0

- Protein context targets 2*K-1 residues for peptide size K, and the `balanced`
  preference is the default ([#219](https://github.com/openvax/isovar/pull/219),
  [details](docs/protein-selection.md)).
