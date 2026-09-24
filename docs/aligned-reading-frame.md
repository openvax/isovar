# Alignment-aware local reading frames

BAM-derived RNA prefixes transfer coding phase through their observed
alignment to the compatible transcript, so an upstream insertion or deletion
is not hidden by trimming RNA and reference strings to equal lengths
([#265](https://github.com/openvax/isovar/issues/265)).

The first retained, mapped RNA base anchors the local frame. Internal query
insertions and exonic reference deletions count as nucleotide differences;
annotated introns do not. Conflicting mappings or transcript-frame hypotheses
yield no translation, not a sequence-only fallback. An included annotated
start codon must remain mapped and intact. This does not infer unobserved
upstream variants, establish a full-length CDS, or distinguish a biological
indel from a sequencing/alignment error.

Manually constructed `AlleleRead` objects without alignment metadata are
matched by sequence alone; supplying alignment metadata is necessary for the
indel-aware guarantee.

Regression coverage includes both strands, repetitive prefixes, and original
Sid CD109 ONT observations with a one-base deletion. RNA translation is checked
separately from a reference-transcript-plus-single-edit prediction. These are
different hypotheses, not interchangeable ground truth.

Coordinates follow the query/reference consumption rules in the
[SAM specification](https://samtools.github.io/hts-specs/SAMv1.pdf).
