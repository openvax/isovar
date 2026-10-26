# COLO829 shared-workflow inventory

This replaces the selected-path computation in the SV3 path/frame report
scripts with the same library workflow used for every retained path and every
resolved event definition. The older `sv3_path_fusion.json` and
`sv3_frame_anchor_audit.json` one-path reports are historical snapshots;
[REPORT.md](REPORT.md), [SV3.md](SV3.md), and `results/` are the new outputs.

## Reproduce

The committed inputs are sufficient for an offline replay. Run from the
repository root, with the version and dependencies recorded in `report.json`:

```sh
python -m isovar sv-audit \
  --manifest reports/hcc1395_colo829/unified_sv_audit/manifest.json \
  --output-dir /tmp/colo829-sv-audit
python -m reports.hcc1395_colo829.build_unified_sv_report \
  --results /tmp/colo829-sv-audit --output /tmp/colo829-sv-audit.md
```

Use `--resume` only for the identical implementation, manifest, input bytes,
indexes and dependencies. Select `--event-id PaulinPaz_2024_v2:COLO829_SOMATIC_SV3_DEL`
for the SV3 subset. The report records a subset selection explicitly.

`build_unified_sv_manifest.py` regenerates inputs from the two original VCFs
and locally installed Ensembl 87. Varcode 10.19.1 owns adjacency interpretation
and all predicted hypotheses; unresolved materialized cDNA remains null.
Every transcript from a gene overlapping a breakpoint +/-10 kb is retained
when its sequence is available. CDS validation failures retain the noncoding
reference sequence and a reason, not an invented frame. The annotation has no
attached chromosome FASTA, so Varcode's unannotated cryptic-exon search is
incomplete; that search is not claimed here. Changing annotation or prediction
versions requires a new manifest and run.

## Scope

The manifest preserves all 54 PaulinPaz 2024 v2 and 134 lifted ValleInclan 2022
VCF records separately, including paired breakends. These are 188 source
records, not 188 independent somatic events. Both RNA orientations are explicit
views; inversions have both adjacencies. Four malformed BND definitions and one
reversed INV interval remain unassessed. Nine symbolic insertion records lack
a resolved oriented adjacency. See [input_problems.json](input_problems.json)
and [#575](https://github.com/openvax/isovar/issues/575).

The SRR8615617 selected-site BAM was extracted after alignment of the entire
115,269,780-pair run. It preserves 500 bases around nominated SV ends and the
small-variant atlas sites; it is not a whole-genome alignment. SV3 also uses the
CNIH3-wide original-read extract and ten original FASTQ ends. Both extracts are
the **same library**, with overlapping fragments. The runner never pools them
or treats them as independent replication. Each source records its sample,
library, culture, role, acquisition limits and checksums.

The original PacBio 2019 tumor and normal RNA sources remain explicit
unavailable-source rows for SV3. On 2026-10-26 the documented `downloads-ap`
endpoint timed out and the alternate `downloads` endpoint returned 404. The
Zenodo/TAGET mirror has polished representatives and counts, not original
molecules. Recovery is tracked in [#574](https://github.com/openvax/isovar/issues/574).
The historical polished cross-check remains separately available in
[colo829_sv3_polished_check.md](../colo829_sv3_polished_check.md); it is not an
original-read source in this manifest.

The workflow records every retained branch and all source limitations.
Discovery budgets can prune paths. Exact matching can miss noisy reads.
A no-candidate result, unavailable normal RNA or bounded search cannot prove
absence. SV3 frame resolution remains in [#563](https://github.com/openvax/isovar/issues/563).
The broader small-variant atlas remains in [COLO829_ATLAS.md](../COLO829_ATLAS.md).

## Outputs

Each event/source/view has a compressed JSON snapshot and SHA-256 receipt.
Snapshots retain reconstruction, all-path assessment, competitive hypotheses,
original witness references and the exact manifest. Exon-alignment records are
stored once in `assessment.exon_alignments`; each path link names its
`exon_alignment_id`. Raw queries, flags, placements and clip quality remain
inspectable without copying the same alignment into every branch.

The protein counts in the report concern annotated-frame candidates, split by
partial versus complete sequence. Exploratory ORFs remain separate. The API
also accepts the existing supplied-fusion input format. Tests exercise the
HCC1395 531-aa complete candidate and partial BCR–ABL1 control, plus synthetic
alternate placements, unresolved branches, reverse-strand CDS, copied records,
competing/coexisting hypotheses and checksum-verified replay.
