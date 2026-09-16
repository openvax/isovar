# Observed coding-fusion hypotheses

These original, unchanged RNA records support coding hypotheses, **not unique
coding transcript assignments or clinically validated antigens**. All available
Ensembl 87 transcript models are retained, including incompatible CDS alternatives.
Inputs, models and original SAM records are pinned by the adjacent manifest.

- **Sid T0, ATP5MG--KMT2A**: Personalis Isofox `Id_14364`, RNA junction
  chr11:118401717 / 118468775. Two distinct, nonduplicate template names
  support the same 105-nt window, Q20 at every base. The observed annotated
  donor start lies 8 nt into the window. The coding hypothesis is
  `MAQFVRNLVEKTPALVNG*`: a mixed junction codon at offset 17, then a stop.
  There are four junction-spanning 8--11-mers. Five compatible transcripts
  have no complete annotated CDS, so the result remains `ambiguous`.
  The donor gene's Ensembl 87 name is ATP5L (current name ATP5MG).
- **K562, BCR--ABL1**: external positive-control RNA, not Sid. Two exact
  sequence-and-mapping windows span BCR exon 14 / ABL1 exon 2. The 240-nt
  sequence supports a 79-aa conditional local translation, including 38
  junction-spanning 8--11-mers. Two donor models give the same frame; a third
  places this join after its CDS. The result remains `ambiguous`, with no
  claimed observed full CDS start. The original pool combines ENCODE Iso-Seq
  products, so per-read library provenance cannot be resolved further.

These are bounded examples, not whole-sample abundance estimates. Peptide
windows are not predictions of binding, tumor specificity or proteome novelty.
An RNA join does not by itself prove a somatic DNA rearrangement.

Sources: [Sid fusion catalogue](https://osteosarc.com/fusions/),
[JBrowse K562 RNA example](https://jbrowse.org/jb2/docs/tutorials/k562_fusions/).
Exact BAM URLs and original records are embedded in each input. No reference
sequence fills an unobserved RNA base and no RNA base or alignment is edited.

Rebuild with `python -m tests.data.fusions.build_coding_examples --sid-bam ...
--k562-bam ... --output tests/data/fusions/coding-corpus`. The regional BAMs
were acquired with indexed `samtools view --no-PG -b -X`: Sid regions
chr11:118401600-118401800 and chr11:118468700-118468900; K562
chr22:23286000-23293000. The public indexes have the BAM URL plus `.bai`.
