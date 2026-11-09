# Intronic insertion model outcomes (#689)

An insertion anchor outside a supplied transcript's exons currently leaks
`ValueError` from `Transcript.spliced_offset`. The shared coordinate helper must
raise its existing `NoExonicReferenceBases` outcome for that case, as it already
does for substitutions/deletions. Reconstruction can then retain the unresolved
model and continue evaluating other supplied ORFs.

Convert only the insertion anchor lookup's `ValueError`, retaining its cause.
Preserve strand-specific insertion offsets, exon-boundary behavior, literal
allele validation and unexpected errors. Do not catch all reconstruction errors
or turn a reference-model failure into zero RNA support.

Use synthetic, sequence-validated PyEnsembl/ORFdata models on both strands to
test intronic and outside-transcript anchors, existing exonic/boundary offsets,
reference-key rejection, and a mixed valid/intronic model assessment with and
without reads, in either model order. Run focused regressions, lint, the full
suite, all CI jobs, then deploy the patch release from clean master.

Primary API reference reviewed: [PyEnsembl Transcript.spliced_offset](https://github.com/openvax/pyensembl/blob/main/pyensembl/transcript.py)
documents an exonic position requirement and `ValueError` otherwise. This fix
changes error classification, not splice or translation biology.
