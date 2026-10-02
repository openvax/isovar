# Anchored flanking-indel comparison (#428)

## Problem and scope

Equal-length clipping before global flank alignment turns a single insertion
into two edits by shifting the distal window boundary. Correct the shared
matcher used by cell attribution and optional partial-read support. Discovery,
fragment allocation, focal-allele quality categories and raw counts stay as-is.

## Comparison contract

- Compare each complete flank outward from the focal allele. Reverse the
  prefix strings so their variant-adjacent ends become the alignment starts.
- Anchor both starts, charge unit substitution/insertion/deletion costs, and
  stop when either sequence ends. The remaining distal overhang is unobserved
  in the other sequence and is free. Never allow a free gap at the focal end.
- Compute this overlap with Edlib prefix alignment in both query/target
  directions. Keep every optimal endpoint pair, including ties; sequence
  length alone does not determine which sequence ends first after indels.
- Use the minimum consumed read length over optimal endpoints for the overlap
  floor, including observed insertions. Use the minimum consumed candidate
  length for `candidate_interval` and `spans_candidate`, so an endpoint tie
  cannot assert certain full-candidate coverage. Export endpoint pairs and
  the compared-read-base count to make these decisions inspectable.
- Check transcript compatibility and A/C/G/T content over the maximum read
  span of tied endpoints, and A/C/G/T content over the maximum candidate span.
  A tie must not hide an incompatible splice path or unknown compared base.
- Keep the focal allele globally aligned and separately anchored. Its paired
  reference competitor uses the same flank costs; ties still reject
  mutant-specific support. Preserve alternatives within the edit tolerance.
- Record the changed comparison policy as `anchored_edit_compatibility.v2`
  in cell and fractional-support exports. Allocation remains v1.

Distal terminal gaps cannot distinguish an indel from incomplete coverage
using sequence alone; they are free overhangs, with conservative span reporting.
This is a local compatibility rule, not an error-probability model.

## Validation and release

Reproduce #428 before the fix. Test left/right internal and anchor-adjacent
insertions/deletions, distal overhangs, partial windows, tied endpoints,
minimum overlap, unknown bases, transcript paths and reference competitors.
Check distances/endpoints against an independent small dynamic-programming
oracle. Exercise actual BAM collection, pooled reconstruction, cell attribution
and fractional scoring on both transcript strands with independently known
proteins. Run existing real Illumina/ONT/PacBio regressions, full lint/tests and
CI; bump to 1.42.1, merge and deploy through the standard script.

## Primary basis

[Edlib's alignment-mode documentation](https://github.com/Martinsos/edlib#alignment-methods)
defines prefix alignment with a fixed start and a free target suffix. Taking
the minimum in both directions gives the either-sequence-end boundary above.
The underlying unit-cost alignment implementation is described in
[Šošić and Šikić, 2017](https://doi.org/10.1093/bioinformatics/btw753).
