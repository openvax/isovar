# ORF start evidence

When [`isovar sv-rna`](sv-rna.md) finds an open reading frame (ORF) crossing an
SV, a key question is where translation would begin. Each exploratory ORF's
`start_evidence` says how well its ATG is supported by annotation and by the RNA,
in four tiers. The JSON and TSV exports carry it unchanged.

A tier describes where a start could plausibly be; it does not show that
translation starts there. Initiation, translation, the RNA's strand and whether
a segment is in a mature transcript remain separate questions. The
classification never changes which candidates exist, and never merges
conflicting isoforms into the strongest tier.

## Tiers

| Priority | Tier | What it means |
| --- | --- | --- |
| 1 | `annotated_CDS_start` | The ATG is the start codon of a supplied same-strand coding transcript, and the ORF reaches the event before stopping. This does not show that translation starts there in the sample. |
| 2 | `five_prime_UTR_ATG` | The ATG sits in the 5′ leader of a supplied same-strand coding transcript. The record says whether, in the reference, it reads as an upstream ORF, an overlapping ORF or an in-frame extension. An ATG alone does not prove an upstream ORF is translated. |
| 3 | `intronic_splice_supported` | The ATG is intronic, and one read covering the whole ORF links it to the event and to qualifying transcript-inclusion evidence (below). Initiation is still not observed. |
| 4 | `sequence_only` | Any other ATG: intronic without that support, antisense, inside a CDS, in a 3′ UTR, noncoding, created by the junction, or in unplaced bases. The subtype is kept. "Outside the supplied transcripts" is not a genome-wide intergenic claim. |

Reading the tiers:

- A lower number means stronger annotation support. Priorities are not
  calibrated probabilities or a recommendation for peptide selection.
- When transcripts or paths disagree on the tier, the priority is null and every
  alternative is kept.
- Read support and read quality are reported separately from the tier.
- Unplaced bases can be soft clips or junction homology, not necessarily
  inserted DNA.
- If the canonical frame stops before the event, it cannot support an
  independent downstream intronic ORF. No splice is chosen just because it
  would remove that stop.

## Priority 3: an intronic start in an observed transcript

An intronic ATG could be translated if the RNA shows its intron segment spliced
into a mature transcript. Isovar checks this on a strand-aware graph of the
supplied exons, observed CIGAR `N` junctions, linked supplementary alignments
and retained segments:

- It distinguishes cryptic donor or acceptor use, exonization, intron retention
  and rearranged intronic segments.
- A genomic deletion or supplementary join alone is not splicing. Retention
  needs both exon/intron boundaries, other annotated splicing on the same read,
  and other reads that splice the intron out.
- The evidence must come from one read that covers the whole ORF (a *full-ORF
  witness*). Its MAPQ and base qualities are checked over the linked interval,
  and missing qualities are reported as unavailable, not low.
- Witnesses are combined across reconstruction modes, not summed.
- Motif-only, wrong-strand, low-quality, unspliced-only and unresolved cases
  never reach priority 3.

The thresholds, failures, witnesses, signal ancestry and library-scoped UMI
support are all exported. The gates and their options are described in
[linking an intronic start](sv-rna.md#link-an-intronic-start-to-an-observed-transcript-path).
Matched-normal comparison and splice prediction on the alternate haplotype are
not performed. They are recorded as `not_assessed` and never stand in for
observed linkage.

`annotate_orf_start` assigns the tiers, and `annotate_orf_inclusion` runs the
priority-3 assessment on retained observations.

## Scientific basis and real-data checks

Calvo et al. distinguish upstream and overlapping ORFs from the main CDS and
show why an upstream AUG does not establish its translation:
[PNAS 2009](https://pmc.ncbi.nlm.nih.gov/articles/PMC2669787/).
For retained-intron candidates, RNA inclusion is a different question from
translation or presentation: [Smart et al., 2018](https://pubmed.ncbi.nlm.nih.gov/30114007/).
Splice prediction is a separate, model-based prior:
[Jaganathan et al., 2019](https://doi.org/10.1016/j.cell.2018.12.015).

Pinned original-read fixtures check three cases: TPST1–CRCP keeps its 30-aa
5′-UTR candidate, PARD3B keeps its intronic and unplaced-start alternatives, and
GABBR1–SLC29A1 keeps its stop-before-join distinction. These are evidence
annotations, not vaccine selection or evidence of antigen presentation.
