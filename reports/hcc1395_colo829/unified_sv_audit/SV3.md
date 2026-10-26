# COLO829 SV3 shared-workflow assessment

1 source-defined events; 8 event/source/orientation runs.

Every retained path is assessed by `isovar.assess_sv_rna_paths` and compared by `isovar.reconcile_sv_rna_hypotheses`. Complete results and original witness provenance are in the linked compressed snapshots. No report-specific frame inference is performed.

Counts below describe retained computational paths, not independent molecules or resolved biological events. The two SRR8615617 subsets are the same library. Opposite orientations and paired VCF records are separate hypotheses, not replication. Exact-sequence compatibility does not validate base quality, establish translation or prove somatic causation.

| Source | Library | Culture | Role |
|---|---|---|---|
| SRR8615617-selected-sites | SRR8615617 | SRR8615617-separately-grown-from-benchmark-DNA | tumor |
| SRR8615617-CNIH3 | SRR8615617 | SRR8615617-separately-grown-from-benchmark-DNA | tumor |
| PacBio2019-COLO829T | PacBio2019-COLO829T | PacBio2019-separate-ATCC-culture-COLO829T | tumor |
| PacBio2019-COLO829BL | PacBio2019-COLO829BL | PacBio2019-separate-ATCC-culture-COLO829BL | normal |

Run outcomes: event_linked_candidates: 4, unassessed: 4.

Missing original PacBio tumor/normal evidence remains open in [#574](https://github.com/openvax/isovar/issues/574). Malformed or sequence-unresolved source definitions remain open in [#575](https://github.com/openvax/isovar/issues/575). An unassessed or bounded no-candidate result is not evidence of absence.

| Event / source / view | Status | Retained paths | Direct event join / compatible join paths | Event-linked partial / complete protein candidates | Limits |
|---|---|---:|---:|---:|---|
| [PaulinPaz_2024_v2:COLO829_SOMATIC_SV3_DEL / SRR8615617-selected-sites / 1-forward](results/f3cdf5f7ee224283800fcfa9.json.gz) | event_linked_candidates | 1 | 1 / 0 | 0 / 0 | selected-site evidence cannot exclude unsampled distant exons or placements, different processed subsets are not independent libraries |
| [PaulinPaz_2024_v2:COLO829_SOMATIC_SV3_DEL / SRR8615617-selected-sites / 1-reverse](results/65d4d8a57fac20b2740ec90e.json.gz) | event_linked_candidates | 1 | 1 / 0 | 0 / 0 | selected-site evidence cannot exclude unsampled distant exons or placements, different processed subsets are not independent libraries |
| [PaulinPaz_2024_v2:COLO829_SOMATIC_SV3_DEL / SRR8615617-CNIH3 / 1-forward](results/8e0cf4fd07b87310d646499b.json.gz) | event_linked_candidates | 80 | 1 / 0 | 0 / 0 | mate_location_unavailable, path_limit, selected-site evidence cannot exclude unsampled distant exons or placements, different processed subsets are not independent libraries |
| [PaulinPaz_2024_v2:COLO829_SOMATIC_SV3_DEL / SRR8615617-CNIH3 / 1-reverse](results/46f75d65f584d95543042590.json.gz) | event_linked_candidates | 80 | 1 / 0 | 0 / 0 | mate_location_unavailable, path_limit, selected-site evidence cannot exclude unsampled distant exons or placements, different processed subsets are not independent libraries |
| [PaulinPaz_2024_v2:COLO829_SOMATIC_SV3_DEL / PacBio2019-COLO829T / 1-forward](results/41e41c33d8a265f944124839.json.gz) | unassessed | 0 | 0 / 0 | 0 / 0 | Original CCS/FLNC RNA unavailable: documented PacBio host timed out 2026-10-26; alternate host 404. Polished representatives are not original molecules. See #574. |
| [PaulinPaz_2024_v2:COLO829_SOMATIC_SV3_DEL / PacBio2019-COLO829T / 1-reverse](results/590b11fdca58e011e62a55b2.json.gz) | unassessed | 0 | 0 / 0 | 0 / 0 | Original CCS/FLNC RNA unavailable: documented PacBio host timed out 2026-10-26; alternate host 404. Polished representatives are not original molecules. See #574. |
| [PaulinPaz_2024_v2:COLO829_SOMATIC_SV3_DEL / PacBio2019-COLO829BL / 1-forward](results/8b05989e8d581f8f7cf98939.json.gz) | unassessed | 0 | 0 / 0 | 0 / 0 | Original CCS/FLNC RNA unavailable: documented PacBio host timed out 2026-10-26; alternate host 404. Polished representatives are not original molecules. See #574. |
| [PaulinPaz_2024_v2:COLO829_SOMATIC_SV3_DEL / PacBio2019-COLO829BL / 1-reverse](results/e5b5e087b3522ce65416df83.json.gz) | unassessed | 0 | 0 / 0 | 0 / 0 | Original CCS/FLNC RNA unavailable: documented PacBio host timed out 2026-10-26; alternate host 404. Polished representatives are not original molecules. See #574. |

Protein candidates retain their annotated frame evidence and completeness. Partial candidates can assume an upstream annotated frame; exploratory ATGs are reported separately and never become annotated protein candidates merely because reads tile them. Multiple local hypotheses can coexist.

Reproduce with `python -m isovar sv-audit --manifest reports/hcc1395_colo829/unified_sv_audit/manifest.json --output-dir NEW_DIRECTORY`, then `python -m reports.hcc1395_colo829.build_unified_sv_report --results NEW_DIRECTORY`.
