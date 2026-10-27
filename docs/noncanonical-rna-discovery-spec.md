# Mutation-supporting RNA beyond canonical coding contexts

## Report review and priority

The released 1.60.1 reports now distinguish local RNA windows, incomplete
proteins, unresolved SV frames and the DDX17 coding-reference exclusion.
Their largest immediate discovery gap is that reference ncORF overlaps are
counted but never reconstructed as mutation-bearing RNA/protein candidates.
This is distinct from lack of translation evidence. The next stage will expand
SV original-read recovery through the existing shared reconciliation workflow.

## First implementation

1. Add a library API that reconstructs mutation-supporting local RNA and assesses
   supplied, validated ORFdata TranscriptORF models, including noncoding parent
   transcripts, without the ordinary annotated-CDS eligibility gate. Reuse
   Isovar assembly, reference matching and witness accounting. Preserve model
   identities, all retained sequence alternatives, no-candidate outcomes and
   explicit reasons; never replace a parent transcript's annotation identity
   with a fabricated canonical CDS.
2. Validate species, assembly, source transcript/path and reference sequence.
   Distinguish an ORF coding-frame hypothesis from observed initiation or
   translation. Keep local/complete sequence support separate from focal-allele
   support. Do not turn reference MS/Ribo-seq observations into mutant evidence.
   Compare RNA candidates against the corresponding reference ORF and preserve
   unresolved nonfocal differences. Ordinary coding models remain competitors.
3. Replay all ALT-supported ncORF overlaps in the COLO829/HCC1395 inventories,
   with each library separate. Pin sources and candidate evidence, add an honest
   discovery funnel, and highlight newly recovered amino-acid-changing candidates
   separately from synonymous contexts and failures. Deliver updated PDFs before
   lengthy release checks, then run lint/tests/CI, merge and deploy.

## Subsequent SV recovery

Use the DDX17 recovery lessons to expand original mates and nearby exon evidence
for the most supported unresolved SV paths. Include soft clips and unplaced ends;
report discovery bounds and competing placements. Feed the same existing shared
SV assessment/reconciliation API rather than adding a report-only protein call.

Primary scientific context: Mudge et al., *Standardized annotation of translated
open reading frames* (2022), doi:10.1038/s41587-022-01369-0; Cao et al., MiSplice
(2020), doi:10.1038/s41467-020-19307-6. RNA and reference ORF annotation nominate
coding candidates; neither establishes mutant translation or presentation.
