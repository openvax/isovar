# COLO829 SV3 protein audit

## Question

Can the 1,213-base original-RNA assembly at the nominated COLO829 SV3
deletion support a specific altered protein sequence?

## Method

1. Reconstruct the event from the pinned original-read BAM, preserving all
   mapping qualities, soft clips, original fragment IDs and assembly limits.
2. Place every assembled ATG-to-stop candidate relative to the RNA junction,
   Ensembl 87 exons/CDS and available ORFdata catalogues. Distinguish direct
   junction support from flanking coverage and full-length phase.
3. Compare candidate amino-acid sequence across the altered join with the
   same start in the unrearranged GRCh38 locus and screen for reference-like
   or ambiguous peptide sequence. An ATG alone is no evidence of initiation.
4. Export reproducible candidate sequences and an evidence table. Call an
   RNA-predicted proteoform only if a mature transcript path, reading frame and
   mutation-specific sequence are supported. Protein evidence is a separate
   expression validation. Otherwise keep the precise unresolved status and
   the sequences as *hypotheses*.
5. Scan the pinned TAGET polished representative FASTQ in both orientations
   for exact junction-spanning nucleotide windows. Treat a reverse-complement
   match to ordinary downstream reference sequence as a decoy, even if that
   isoform has tumor full-length counts.

## Acceptance

- The targeted SV3 audit is runnable on a supplied indexed BAM with pinned
  reference and annotation, and a regression uses the small COLO829 fixture.
- Multiple reading frames remain separate. No whole-context read is required,
  and the 533 flanking voting fragments are not promoted to SV-isoform support.
- FASTA/TSV/JSON identifiers link each candidate to the original reads,
  junction, start basis, reference comparison and evidence limits.
- The report states whether this event yields a supported proteoform, an
  exploratory ORF, or no call, with a test for that decision.
