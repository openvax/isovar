# Bind a whole-observation seed to its source placement (#629)

The variant-free assembler can join two alternative placements of one original
read because the seed is not yet in the extension's used-observation choices.
Carry an optional actual seed observation through the shared bidirectional
extension function, including its reciprocal orientation. The whole-observation
caller supplies it; junction-core SV callers retain their existing behavior.

Do not combine alternative placements to extend a read with itself. Keep both
separate hypotheses under the existing path budget. A genuinely independent
overlapping read must still extend the sequence. Add both-strand/order and
independent-read regressions, and run existing variant-free/SV suites, full
tests and exact-head CI before merge and clean-master deployment. No new
dependency, CLI flag or alternative reconstruction algorithm is needed.
