# Changelog

Behavior changes that can alter results or break callers, by release. Patch
releases that only fix bugs or add fixtures are omitted; see the
[commit history](https://github.com/openvax/isovar/commits/master) for those.

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
