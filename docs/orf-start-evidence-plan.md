# ORF start evidence and splice inference (#350)

## First PR: start-origin annotation

Add a public, transcript-specific annotation of each event-crossing ATG ORF.
Carry it from RNA reconstruction into JSON and TSV without changing candidate
identity or combining conflicting isoforms into the strongest tier. Rank a
same-strand 5′-UTR ATG above an unsupported intronic ATG, while keeping
initiation, translation, RNA polarity, and transcript inclusion separate.
Preserve the existing reference-ORF comparisons and stop-before-junction gate.
Store the synthetic fixture builders alongside tests and exercise reconstruction
through export with assembly enabled and disabled.

The first PR does **not** infer splice paths or promote an intronic ATG on the
basis of coverage, motifs, or an unrelated splice junction. Tier 3 is reserved
until the linked-evidence assessment below is implemented. Issue #350 remains
open for that work.

## Tiers

| Priority | Tier | Evidence and limits |
| --- | --- | --- |
| 1 | `annotated_CDS_start` | This ATG occupies the exact start codon in a supplied same-strand coding transcript. The reported ORF must itself reach the event before stopping. This does not establish initiation in the sample. |
| 2 | `five_prime_UTR_ATG` | The ATG occupies consecutive spliced positions in the leader of a supplied same-strand coding transcript. Record whether the reference reading is an upstream ORF, an overlapping ORF, or an in-frame extension; ATG presence alone is not proof of a uORF. |
| 3 | `intronic_splice_supported` (reserved) | Intronic ATG linked to qualified transcript-inclusion evidence on the same RNA path; initiation remains unobserved. Never emitted by the first PR. |
| 4 | `sequence_only` | Unsupported intronic, antisense, internal-CDS, 3′-UTR, noncoding, junction-created, or unplaced ATG; retain the subtype. Outside supplied transcripts is not a genome-wide intergenic claim. |

Priorities are an annotation prior, not a calibrated probability or a peptide
selection recommendation. A disagreement between transcript or path assessments
has null priority and retains all alternatives. Missing annotation in an older
result stays unavailable. Read support and quality remain independent axes.
Unplaced bases can be clips or junction homology, not necessarily DNA insertion.

## Follow-up: observed splice linkage

1. Build a strand-aware graph from supplied exons, CIGAR N junctions, linked
   supplementary alignments, retained segments, and DNA breakends. Keep competing
   paths and distinguish alignment orientation from biological RNA polarity.
2. Classify cryptic donor/acceptor use, exonization, exon skipping, intron
   retention, fusion/rearranged exons, and unresolved paths. A genomic deletion
   or supplementary join alone is not proof of splicing. Retention needs both
   boundaries, processed-transcript context, and competing spliced support;
   coverage alone can be pre-mRNA.
3. Require a full-interval witness linking the ATG-containing segment, splice
   mechanism, and event. Assess MAPQ and base qualities over that linked interval;
   missing qualities are unknown. Export thresholds, failures, read identities,
   sample/source, templates, signal ancestry and library-scoped UMI support.
   Union witnesses across reconstruction modes instead of adding their counts.
4. Compare original RNA and competing isoforms, including matched normal where
   available. Absence of informative coverage is unknown, not absent expression.
5. Treat splice motifs and model scores as predictions on the **actual alternate
   haplotype**, with flanks, model/version, coordinates, and limitations recorded.
   They must never stand in for observed linkage or assume calibration on SVs.
6. Transfer annotated frames first; retain independent ATG hypotheses separately.
   Report the first altered amino acid and first stop. A canonical frame that
   stops before the event cannot support an independent downstream intronic ORF;
   never choose a splice just because it removes that stop.

Acceptance tests must differ only in observed linkage and assert a tier change.
Motif-only, wrong-strand, low-quality, unspliced-only and unresolved cases must
not enter tier 3. Add paired tests for assembly enabled/disabled and preserve
allele alternatives in PARD3B rather than pooling their support.

## Scientific basis and real-data checks

Calvo et al. distinguish upstream and overlapping ORFs from the main CDS and
show why an upstream AUG is not sufficient to establish its translation:
[PNAS 2009](https://pmc.ncbi.nlm.nih.gov/articles/PMC2669787/).
For retained-intron candidates, RNA inclusion is a different question from
translation or presentation: [Smart et al., 2018](https://pubmed.ncbi.nlm.nih.gov/30114007/).
Splice prediction is a separate model-based prior:
[Jaganathan et al., 2019](https://doi.org/10.1016/j.cell.2018.12.015).

Pinned original-read fixtures should retain TPST1–CRCP's 30-aa 5′-UTR candidate,
PARD3B's intronic and unplaced-start alternatives, and the GABBR1–SLC29A1
stop-before-join distinction. These are evidence annotations, not vaccine
selection or evidence of antigen presentation.
