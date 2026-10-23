# Original-read SV junction audit

## Question

The broad SV screens count alignments within 500 bp of DNA breakends. Which
candidate links actually resolve an RNA junction, and how many original
fragments support the junction versus a complete observed path?

## Evidence rules

1. Keep source DNA record IDs and screen counts separate. Repeated breakend
   records and two benchmark truth sets may describe the same RNA witness.
2. For an oriented adjacency, retain every usable original RNA alignment
   regardless of MAPQ, including observed soft-clipped bases and actual
   primary/supplementary records. Use dense support to count distinct
   fragments crossing each nominated RNA join. Report MAPQ >= 20 as a
   sensitivity subset, not a validation gate: a clipped partner of a real
   chimeric read can map ambiguously to the unrearranged genome. A mate pair
   is context only and never becomes direct junction support. This SV-specific
   collection choice does not change the ordinary SNV mapping-quality filter.
3. Freeze a small original-read fixture for the three HCC1395 events that
   show an orientation-compatible join. The existing COLO829 selected-site BAM remains
   the paired-read fixture.
4. Extend supported seeds through other observed RNA fragments and retain
   assembly limits. No single read must span the complete assembled context.
   Regional extension is an RNA context hypothesis; count local branch and
   junction support separately from reads that cover other parts of the path.
   Keep `references=[]` in this read-geometry audit. It cannot test Isovar's
   reference-anchored frame transfer or call an SV-derived protein. Ordinary
   splice alternatives and DNA-to-RNA causation remain unresolved.

## Next method: exon-anchor-first reconstruction

Tracked in [Isovar #555](https://github.com/openvax/isovar/issues/555).

The present engine seeds at a placed junction or a clip precisely at the DNA
breakpoint and extends through reads that share a placed genomic base. That
misses mature RNA whose DNA breakpoints lie in excised introns, and it cannot
bridge two clipped read clouds joined only by sequence. A general SV search
must start with an aligned segment inside a reference exon, retain the full
original query and low-MAPQ/soft-clipped partner, and build a sequence graph
from overlaps without requiring the unexpected tail to map to a nominated
DNA breakpoint. The aligned exon gives direction and a local origin; the
unplaced bases are data, not missing reference sequence.

Candidate graph paths should be aligned against rearranged and unrearranged
genome/transcript models and paired RNA from the normal when available.
Exon skipping, cryptic exons, intron retention, inserted sequence, tandem
duplication, inversions, translocations, read-through and ordinary splicing
must compete. Remap original read pairs to each path and count distinct
fragments spanning a *path-discriminating* junction separately from mate
context, compatible exonic coverage and unphased extension reads. Report
multiple placements and unresolved branch phase rather than choosing one
ORF. Only after the RNA path and frame are resolved should it produce an
SV-derived mutant protein. This requires a separate implementation and
normal/reference competition; the selected-event audit below is the
read-geometry baseline, not that complete method.

## Deliverables

- A reproducible script and TSV of audited candidate events, with explicit
  BAM identity and limits.
- A regression test against source-derived original reads, including the
  low-MAPQ COLO829 split-junction fragments that the MAPQ >= 20 sensitivity
  subset misses.
- Report text and PDFs showing direct-junction counts and bounded assembly
  context, with MAPQ retained as alignment provenance.
- [Isovar #553](https://github.com/openvax/isovar/issues/553) for an
  intermediate, path-specific local-context metric;
  entire-path exact support is too strict for noisy long reads, while shared
  junction support is too broad to assign a transcript isoform.
