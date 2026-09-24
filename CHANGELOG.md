# Changelog

Behavior changes that can alter results or break callers, by release. Patch
releases that only fix bugs or add fixtures are omitted; see the
[commit history](https://github.com/openvax/isovar/commits/master) for those.

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
