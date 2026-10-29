# Validated SV source interpretation (#575, #588)

Preserve source bytes, checksums, raw rows and caller/PCR metadata. Put source
interpretation in a typed library module; keep Varcode responsible for VCF
breakend geometry. One API and CLI should inspect standard and explicitly named
source profiles, returning resolved, unresolved and invalid records separately.
The existing manifest preparation must consume that interpretation instead of
maintaining a second parser.

The COLO829 Paulin Paz v2 file stores four literal breakend alleles in END while
ALT is symbolic BND. A checksum-pinned source profile can interpret this exact
encoding, retaining the original fields and a named transformation. Validate
partner contigs/orientations through Varcode and compare with independently
published Minda amended rows where they overlap. Do not match amended event IDs:
their numbering differs. Do not repair the zero-width inversion from SVLEN or
invent sequence for symbolic insertions.

The HCC1395 PCR-labelled file contains 1,118 records and 32 true PCR flags. The
author's GitHub ensemble contains 1,096 rows; all shared core columns agree, but
939 INFO fields differ and 22 rows were added. The GitHub intervals lack END;
the PCR-labelled file inserts END=POS. Caller IDs can contain competing endpoint
claims (Minda_192's endpoint difference is seven bases from SVLEN), not a validated
consensus coordinate. Preserve those claims without turning them into resolved
RNA inputs. Feed only defensible oriented adjacencies to the existing audit.

Validate all actual source records, representative malformed and ambiguous
inputs, PCR tri-state semantics, changed checksums, and CLI/API parity. Publish a
compact source inventory with exact provenance and remaining unresolved reasons.
Run the shared RNA workflow for newly defensible events with its existing source
separation and budgets. Full lint/tests, self-review, green CI and clean-master
PyPI deployment remain required.

Primary references: VCF 4.5 specification (samtools/hts-specs, 25 February 2026),
SVScope paper DOI 10.1186/s13059-026-04076-0, Zenodo 18688406 and 13917379,
Goatofmountain/SVScope revision 0ad36860326490e94b8567276aa54535a79e5b26,
KolmogorovLab/minda revision 6e13ee5120280d7a7fb6deef8c773a2d4cc41c8b.

Preparation review: centralize the transcript-model selection, oriented view
construction and Varcode hypotheses in `prepare_sv_source_record`. The COLO829
and HCC1395 manifest scripts should only pin files, select libraries/budgets and
serialize this result. Preserve source contig names; resolve annotation names
once per context. Missing annotation contigs leave a limited RNA nomination,
not an invented intergenic effect or a dropped viral/decoy event.
