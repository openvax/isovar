# SV RNA and protein reconstruction plan

## Goal

Produce RNA supported SV protein hypotheses from original HCC1395 and COLO829
reads. A protein needs an observed RNA path, a versioned coding or ncORF frame
anchor, and local evidence for the path branch. Keep DNA event linkage, RNA
sequence support, frame, and measured protein expression separate.

## Implementation

1. Extend candidate paths through exact sequence overlap when a read's placed
   bases cannot overlap the current path end. Start from a placed exon anchor;
   retain clipped and inserted bases, strand, fragment identity, and competing
   placements. Bound the graph with existing read and path limits. Never add a
   reference base to an RNA path.
2. Count support over local windows containing the novel join and branch
   distinguishing decisions. Report shared direct joins and branch specific
   observations separately; preserve full-path support as a separate measure.
3. Supply pinned transcript models for source replays, transfer frame only
   from exact coding anchors, and export candidate protein sequences with
   observed support and unresolved alternatives.

The current implementation covers sequence-only exon clips, both fusion
orientations, a 30-base unplaced insertion, ordinary splice and repeated
placement controls. The HCC1395 EIF3K–CYP39A1 source replay now yields a
complete 531-amino-acid coding candidate from original short reads, with an
independent ONT junction check. COLO829 SV3 still has no resolved coding path.
Tandem duplications can stop at `repeated_genomic_position`, and path-specific
support remains tracked by #553; #555 stays open for those cases.

## Verification

Synthetic controls cover an intronic fusion joined by clipped read clouds,
unplaced insertion, ordinary splice competitor, repeated overlap, and low MAPQ
partner. Replay pinned original HCC1395 and COLO829 records. Run lint, tests,
CI, then merge and deploy each release PR. Do not promote a candidate as a
resolved proteoform when transcript orientation, frame, or branch support is
unresolved.
