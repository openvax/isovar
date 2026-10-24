# SV3 RNA path and frame evidence (#563)

## Aim

Reassess the COLO829 SV3 candidate using the pinned original SRR8615617 reads, including the full-source target-screen fixture. Combine direct junction, clipped-sequence, overlapping-read, mate, and versioned exon evidence without counting one SRA spot twice or turning a mate link into sequence continuity. Determine whether any observed path supports a CNIH3 reading frame or distinguishes the existing 37-aa and 162-aa exploratory ORFs.

## Method

1. Replay the source candidate BAM with dense support and exon-anchored clipped bases. Retain all bounded joined paths and record assembly limits.
2. For each exact SV3 junction-window spot, inspect both original read ends, mapping alternatives, soft clips, and overlap with the retained joined path. Compute separate direct-sequence, same-end clipped, mate-location, and exact exon-anchor evidence; never sum unlike evidence into a single coverage count.
3. Compare each exploratory ORF to the observed nucleotide paths. Report start/frame-discriminating evidence separately from RNA evidence shared by all translations of the same path. Only a continuous sequence path to a versioned exon and a justified frame may become an event-linked protein candidate.
4. Publish a reproducible machine-readable table and short interpretation with source IDs, parameters, boundedness, and limitations. Add focused regression tests using the pinned small source fixture.

## Acceptance

The report reruns locally without downloading the full SRA; correctly deduplicates spots across mates and secondary/supplementary records; never treats a mate location as a sequenced bridge; and states clearly whether the evidence resolves a mature RNA path and frame. The complete Isovar lint and test scripts pass. If the original reads leave the frame unresolved, the PR advances #563 but does not close it.
