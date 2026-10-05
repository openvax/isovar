# Preserve distinct original SV alignment records

The focused Sid TPST1::CRCP bulk subset contains 32 collisions under the old
(QNAME, FLAG, RNAME, POS, CIGAR, RG) key. The first pair differs in mate position,
template length and alignment tags. Default collection silently drops one;
dense collection rejects the input. This reproduces #438 and blocks #433.

Hash the complete original SAM alignment text for both indexed deduplication and
serialized record IDs. Exact repeat fetches still collapse, while distinct
payloads remain available to the existing segment, placement and cell/UMI
ambiguity logic. Keep contig-header-order invariance and source-scoped fragment
counts. Do not reinterpret alternate records as independent fragments.

Add regressions for differing mandatory/optional payload fields, identical
repeat records, ambiguous labels and both default/dense reconstruction. Check the
real bulk subset independently, then lint, full tests, CI, merge and release a
patch before replaying the focused audit under the corrected engine.

The alignment payload follows the [SAM specification](https://samtools.github.io/hts-specs/SAMv1.pdf).
Record/observation IDs change once; fragment and segment identities do not.
