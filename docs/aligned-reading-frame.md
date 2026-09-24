# Reading frames from aligned reads

To translate an RNA sequence, Isovar has to know which bases begin codons. It
takes the frame from the annotated transcripts the reads are compatible with,
and carries it along each read's own alignment. It does not line up the RNA and
reference as plain strings.

This matters when the RNA has an indel upstream of the variant. A one-base
deletion shifts the frame of everything after it. Trimming the RNA and reference
to equal lengths and comparing them would hide that shift
([#265](https://github.com/openvax/isovar/issues/265)).

## How the frame is carried

- The first retained, aligned RNA base anchors the frame.
- Each inserted or deleted base in an exon counts as one difference toward
  `max_transcript_mismatches` and moves the frame accordingly. Annotated
  introns (CIGAR `N`) are skipped without counting.
- If a read's alignments, or the frames of its compatible transcripts,
  disagree, there is no translation. Isovar does not fall back to matching
  sequence alone.
- If the sequence includes an annotated start codon, that codon must still be
  aligned and intact.

Reads you build by hand as `AlleleRead` objects without alignment metadata are
matched by sequence alone. Supply alignment metadata to get the indel-aware
behavior.

## What this does not do

It does not infer variants upstream of what the reads show, establish a
full-length coding sequence, or tell a biological indel from a sequencing or
alignment error.

## Tests

Regression tests cover both strands, repetitive prefixes, and original Sid
CD109 ONT reads with a one-base deletion. The translated RNA is checked
separately from a prediction made from the reference transcript plus a single
edit. These are different hypotheses, and neither is ground truth for the other.

Coordinates follow the query and reference consumption rules of the
[SAM specification](https://samtools.github.io/hts-specs/SAMv1.pdf).
