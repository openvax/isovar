# Shared observed-RNA frame assessment (#583)

Move transcript coordinate models and CDS-anchor translation out of the large
SV reconstruction module. Keep existing SV JSON and scientific behavior stable;
SV discovery chooses event-related departures, while the common frame assessor
operates only on an observed nucleotide path and explicitly supplied models.
Use named, typed immutable frame results internally, with small conversion
helpers for the existing JSON schema. Preserve unresolved coding and noncoding
alternatives and sequence/placement distinctions.

For intronic literal discovery, add model-relative retention boundaries: an
observed contiguous genomic path can leave an annotated coding exon into its
intron without a CIGAR N gap. Transfer frame only downstream of an exact coding
anchor. Never fill a missing RNA base from reference. Annotated frames remain
hypotheses, not evidence of initiation or translation. A focal-allele interval
must overlap the translated path to count as an allele-containing protein;
upstream stops are explicit outcomes. Keep partial frame candidates and ordinary
transcript competitors; synonymous/reference-identical contexts remain distinct
from changed proteins.

This PR is the reusable frame layer for the next original-read recovery and
allele-to-path consumer. It does not claim to complete #583 or reclassify the
old DDX17 tiling as one phased transcript. The following consumer must use the
same original-RNA observation geometry, technical-sequence policy and bounded
extension code as SV reconstruction, with ALT and REF seed adapters and scoped
witness counts.

Self-review and validation: unchanged SV snapshot/translation tests; synthetic
plus/minus retention and cryptic paths; true frameshift versus upstream stop;
alternative/noncoding models; no coding anchor; partial and complete starts;
no reference-filled sequence. Run lint, full tests and exact-head CI before
merging, then the clean-master deploy gate.

Primary scientific references checked 2026-10-29:
- MiSplice: https://www.nature.com/articles/s41467-020-19307-6 — distinguishes
  new splice junctions and allele-associated splicing.
- Smart et al.: https://doi.org/10.1038/nbt.4239 — retained introns can yield
  peptides, but coverage alone is not sample-specific translation evidence.
- SAM specifications: https://samtools.github.io/hts-specs/SAMv1.pdf and
  https://samtools.github.io/hts-specs/SAMtags.pdf — alignment coordinates and
  CIGAR skips must retain their original reference context.
