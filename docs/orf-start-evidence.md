# ORF start evidence

Each event-crossing ATG ORF from `isovar sv-rna` carries a transcript-specific
start-origin annotation (`start_evidence`), produced by `annotate_orf_start` and,
for intronic starts, `annotate_orf_inclusion`. It is carried into the JSON and TSV
exports without changing candidate identity or merging conflicting isoforms into
the strongest tier. Initiation, translation, RNA polarity and transcript inclusion
remain separate questions. A canonical frame that stops before the event cannot
support an independent downstream intronic ORF, and no splice is chosen just
because it removes that stop.

## Tiers

| Priority | Tier | Evidence and limits |
| --- | --- | --- |
| 1 | `annotated_CDS_start` | This ATG occupies the exact start codon in a supplied same-strand coding transcript. The reported ORF must itself reach the event before stopping. This does not establish initiation in the sample. |
| 2 | `five_prime_UTR_ATG` | The ATG occupies consecutive spliced positions in the leader of a supplied same-strand coding transcript. Record whether the reference reading is an upstream ORF, an overlapping ORF, or an in-frame extension; ATG presence alone is not proof of a uORF. |
| 3 | `intronic_splice_supported` | Intronic ATG linked by one full-ORF witness to the event and to qualified transcript-inclusion evidence on the same RNA observation; initiation remains unobserved. See [splice-linked inclusion](sv-rna.md#link-an-intronic-start-to-an-observed-transcript-path). |
| 4 | `sequence_only` | Unsupported intronic, antisense, internal-CDS, 3′-UTR, noncoding, junction-created, or unplaced ATG; retain the subtype. Outside supplied transcripts is not a genome-wide intergenic claim. |

Priorities are an annotation prior, not a calibrated probability or a peptide
selection recommendation. A disagreement between transcript or path assessments
has null priority and retains all alternatives. Missing annotation in an older
result stays unavailable. Read support and quality remain independent axes.
Unplaced bases can be clips or junction homology, not necessarily DNA insertion.

## Splice-linked inclusion (priority 3)

An intronic start is assessed against a strand-aware graph of supplied exons,
CIGAR N junctions, linked supplementary alignments and retained segments. Cryptic
donor/acceptor use, exonization, intron retention and rearranged intronic segments
are distinguished; a genomic deletion or supplementary join alone is not proof of
splicing, and retention needs both boundaries, processed-transcript context and
competing spliced support. The linked interval must come from one full-ORF
witness, with MAPQ and base qualities assessed over that interval; missing
qualities are unknown, not low. Thresholds, failures, witnesses, signal ancestry
and library-scoped UMI support are exported, and witnesses are unioned across
reconstruction modes rather than summed. Motif-only, wrong-strand, low-quality,
unspliced-only and unresolved cases never reach priority 3.

Matched-normal comparison and splice prediction on the alternate haplotype are
not performed; they are recorded as `not_assessed` and never stand in for
observed linkage.

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
