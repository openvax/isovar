# Alignment-aware local reading frames

Starting in 1.16.0, BAM-derived RNA prefixes transfer coding phase through
their observed alignment to the compatible transcript. An upstream insertion
or deletion is no longer hidden by trimming RNA and reference strings to
equal lengths (#265). This can change candidate proteins and their support.
Splice-path classification and assembly defaults are unchanged.

The first retained, mapped RNA base anchors the local frame. Internal query
insertions and exonic reference deletions count as nucleotide differences;
annotated introns do not. Conflicting mappings or transcript-frame hypotheses
yield no translation, not a sequence-only fallback. An included annotated
start codon must remain mapped and intact. This does not infer unobserved
upstream variants, establish a full-length CDS, or distinguish a biological
indel from a sequencing/alignment error.

Manually constructed `AlleleRead` objects without alignment metadata retain
the historical sequence-only API behavior. Supplying alignment metadata is
necessary for the indel-aware guarantee; public signatures are unchanged.

Regression coverage includes both strands, repetitive prefixes, and original
Sid CD109 ONT observations with a one-base deletion. RNA translation is checked
separately from a reference-transcript-plus-single-edit prediction. These are
different hypotheses, not interchangeable ground truth.

The 49-case pinned expansion corpus retains all top protein sequences and
outcome labels. Six cases change intermediate translation counts (EXOC4,
H1-2, NME1, NTF3, VPS72, ZNF764). NME1's unchanged top sequence has 40 rather
than 42 supporting observations: two indel-bearing observations no longer
translate into that protein. Input BAMs and reference assets are unchanged.

Coordinates follow the query/reference consumption rules in the
[SAM specification](https://samtools.github.io/hts-specs/SAMv1.pdf).
