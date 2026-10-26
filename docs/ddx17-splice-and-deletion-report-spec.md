# DDX17 RNA-path investigation and deletion-aware report alignment

1. Reproduce why COLO829 DDX17 chr22:38489856 GA>G has ALT RNA but no
   local protein. Inspect every retained original allele witness, its mate
   and clipping/splice evidence; compare against pinned coding and noncoding
   transcript models. Recover broader original sequence when available and
   distinguish an untested splice hypothesis from a negative observation.
2. Correct reference/RNA comparisons for supported in-frame deletions.
   Anchor the gap to the transcript-specific variant effect, align the
   unchanged suffix, keep reference residues neutral, and mark the deleted
   residue interval with an orange vertical boundary. Preserve ambiguous
   comparisons instead of attributing unexplained differences.
3. Produce and visually inspect updated PDFs and deliver them before lengthy
   release checks. Add regression coverage for the deletion comparison and
   traceable DDX17 evidence. Record remaining limitations in issues, then
   validate, merge, and release through the normal repository gates.
