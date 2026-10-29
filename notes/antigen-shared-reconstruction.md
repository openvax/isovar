# Shared discovery and reference resolution

Review before the catalogue-scale survey reproduced two remaining isolated paths:
viral discovery discards independent parents sharing placement geometry (#643),
and antigen reference construction retains a private contig alias implementation
that silently resolves ambiguity (#644).

Use the existing generic observation/discovery engine for each viral contig,
preserving per-contig budgets, original parents and full-input support. Count
unassemblable target segments separately so they do not spend the observation
budget. Use one Varcode resolver per antigen reference context for exact names,
aliases and mitochondrial coding classification; retain structured unresolved
diagnostics. Tests cover independent versus alternative placements, rejected
observations, ambiguity and exact viral/decoy/accession names. No new biological
model or coding annotation is introduced.
