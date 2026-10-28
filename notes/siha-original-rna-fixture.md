# Original SiHa RNA fixture (#531)

Retain a bounded original-read selection from the already acquired public
SRR12650977 competitive alignment. Select complete query-name groups whose
recorded placements are all on HPV16 and support the previously reported local
path. Preserve every selected SAM record, mate, flag, quality and alignment tag.
Remove only unused header sequence entries; retain the complete source header
and exact source BAM/result/reference hashes separately.

Use the existing sequence-pinned HPV16 FASTA for a small custom-reference replay.
Do not claim that this slice repeats whole-library competition: retain the
separate historical full-input receipt (26 paths, 71 target fragments) and its
human-competitive placements. Pin the ORFdata catalogue and source CDS identities.
Assert the existing local RNA/coding prediction, bounded support and ambiguity;
do not infer protein expression, a full viral transcript or viral absence.

Provide a deterministic extraction recipe, hashes, public-source attribution and
redistribution policy, plus one offline API/CLI regression. No reconstruction
algorithm or new biological sequence engineering is part of this change.
