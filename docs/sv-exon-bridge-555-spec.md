# Exon anchored SV path completion (#555)

## Goal

Turn an observed sequence bridge into a mature RNA and frame call only when original reads connect the bridge to a versioned reference exon. Preserve the direct join count, exact local branch count, exon link, and frame evidence as distinct observations.

## Work

1. Replay selected HCC1395 and COLO829 original reads with versioned exon models and exon clip discovery in both bounded and dense modes. Keep all MAPQs and report placement uncertainty rather than suppressing low MAPQ join sequence.
2. Test the overlap path through cryptic or inserted sequence, skipped exons, repeated anchors, inversion, tandem duplication, ordinary splice and readthrough. Ensure a path is assembled from observed read bases only; a full length single read is not required.
3. For SV3, examine original RNA beyond the 500-base breakend windows for a CNIH3 exon to join sequence link. Require an exact reference exon anchor and an observed RNA path from that anchor through the SV join before transferring a reading frame. If the original material does not provide it, leave protein unresolved.
4. Regenerate the source table with separate local join, mature exon path, frame, and protein columns. Keep the candidate FASTA and audit artifacts reproducible.

## Decision rules

An intronic genomic location alone cannot supply a coding frame. An ATG within a partially assembled path is an exploratory ORF, not evidence of initiation. A low MAPQ alignment can support read sequence while leaving genomic placement ambiguous. Zero single-read matches to a 1,213-base path does not contradict a shorter observed join.
