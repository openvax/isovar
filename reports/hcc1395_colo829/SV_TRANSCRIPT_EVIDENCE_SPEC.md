# Original RNA support for SV transcript hypotheses

## Question

For the pinned HCC1395 and COLO829 structural variant candidates, report
which RNA joins and transcript paths can be reconstructed, how many original
RNA fragments support each claim, and how their sequence differs from the
corresponding unrearranged reference. Explain what prevents an SV3 proteoform
call.

## Evidence rules

Count original fragments once per junction, independent of how many assembled
paths reuse that junction. Report MAPQ as placement sensitivity rather than a
required junction quality gate. A local path count requires one original RNA
observation to match an exact sequence and compatible placements across the
specified window; mates and flanking voters do not establish the path. State
whether the interval distinguishes the retained competing paths. Keep partial
RNA paths, mature transcript models, frame predictions, and protein evidence
as separate statuses. Do not infer a reference protein from an intronic
genomic translation.

## Deliverables

1. A reproducible event table from the pinned original-read fixtures, with
   direct junction counts, local path counts and limits, MAPQ sensitivity,
   observed junction sequence, and explicit reference comparison.
2. A human-readable assessment including the larger independent HCC1395
   EIF3K--CYP39A1 fusion, with source provenance kept separate.
3. Regression tests for path-specific support with competing isoforms and for
   the pinned HCC/COLO candidates. Run lint, tests and CI before release.
