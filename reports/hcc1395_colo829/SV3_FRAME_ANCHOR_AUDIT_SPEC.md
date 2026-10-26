# SV3 regional frame-anchor evidence (#563)

## Question

Do any original COLO829 RNA reads connect the exact SV3 junction to a
versioned CNIH3 exon or coding frame? Report other exon-to-region links as
alternative paths, not as support for the SV3-specific RNA product.

## Inputs and method

- Reuse the pinned `SRR8615617` CNIH3 candidate BAM, Ensembl 87 models, and
  `sv3_path_fusion.json` joined sequence. The BAM is a full-source target-screen
  subset, not every RNA alignment in the original SRA.
- Scan every candidate alignment without a MAPQ gate. Match original query
  sequence in either orientation to the joined path for at least 20 exact bases.
  Deduplicate placements by SRA spot and mate; preserve raw alignment details.
- Identify exon placement using at least 15 aligned bases in a versioned exon.
  Distinguish a same-end clipped sequence link from a mate link and from a
  junction-spanning read. Record clipped-base quality and whether that spot
  contains the exact 43-base SV3 join.
- Locate links relative to the SV3 join and reference blocks. Show which
  exons have an annotated CDS and where its start lies. Compare with the
  existing junction witnesses, ORF starts, and polished-isoform screen.

## Interpretation and verification

An exon read matching a downstream reference block can support a competing
splice path without placing the exon on the SV3-specific RNA molecule. Neither
an ordinary reference overlap, mate position, nor an exploratory ATG transfers
a frame. Report search limits, exact read IDs, and unresolved initiation and
translation. Regenerate JSON and Markdown offline; test the source-derived
links and a small synthetic case where exon and junction evidence stay distinct.
