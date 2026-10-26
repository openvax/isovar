# Sid neoORF report: DNA event → RNA transcript → protein sequence → support

The strongest leads in the focused screen are the altered **TPST1–CRCP upstream
ORF**, **SPAG1, TECPR1 and PIP5K1A frameshift windows**, and exploratory
**FOXO3–STRADA/CCDC47 ORFs**. KTN1 has a shorter shifted tail with a stop.
The newly recovered **FTH1–FTL** sequence is a conventional in-frame fusion
ORF lead with an annotated FTH1 start and the normal FTL stop; its exact
complete-ORF support is two ONT T2 fragments. Its coding anatomy is clearer
than the aggregate ambiguous-model label alone suggests.
These are RNA sequence observations and coding hypotheses; none has demonstrated
protein translation or antigen presentation in this analysis.

This report includes **every entry in the selected SV/indel/splice search**,
including unknown DNA origins, unresolved transcripts and zero-candidate
outcomes. It covers **23 SV geometries and 32 small-variant nominations**,
not all possible neoORFs in Sid's tumors. The
[55 event sheets](data/sid-neoorf-event-report/README.md) contain original
alleles/calls, DNA support, transcript anatomy, RNA outcomes, every protein
sequence/window and support by product. Their index includes the unresolved
entries as well as the leads.

The dense recovery contributes 44 completed views, including both ONT IGKC
orientations completed in [#470](https://github.com/openvax/isovar/issues/470).
Their original timeout and interrupted outcome remain preserved as historical
failures with unknown counts. The full 100-view selected SV outcome ledger is
retained; source-specific support is not pooled across products or orientations.

## Complete, downloadable sequences and evidence

- [Every protein/window, FASTA](data/sid-neoorf-event-report/proteins.fasta):
  1,704 event-specific sequence groups, including 723 in-frame control groups.
- [SV ORF nucleotide sequences, FASTA](data/sid-neoorf-event-report/sv-orfs.fna):
  867 groups with exact reconstructed ORF nucleotides.
- [Support table, TSV](data/sid-neoorf-event-report/support.tsv): all 1,785
  source-specific hypotheses: 895 SV ORFs, 116 frameshift windows and 774
  in-frame windows.
- [Full event/placement ledger](https://github.com/openvax/isovar/releases/download/v1.54.0/sid-neoorf-events-1.54.0.json.gz):
  original DNA provenance, transcript models, per-product outcomes, all
  6,165 SV ORF occurrences with genomic RNA blocks and junctions,
  sequence identities, stop/start/frame flags and input hashes.
- [File checksums](data/sid-neoorf-event-report/files.json).

The complete JSON ledger is a versioned release attachment. Its
[download record](data/sid-neoorf-event-report/downloads.json) pins the unchanged
bytes, size and original report manifest. The portable report archive retains
the complete original report with local links.
The matching release also provides a [download verification receipt](https://github.com/openvax/isovar/releases/download/v1.54.0/completed-report-release-download-verified.json)
for the published ledger.

Alternative starts, overlapping windows, synonymous RNA sequences and the two
PARD3B geometries reuse evidence. These totals count hypotheses, not independent
neoORFs. The protein FASTA gives the **entire reported ORF or coding window**;
local indel windows do not reconstruct an entire full-length gene protein.
Their full altered mature transcript and full protein sequence remain unknown.

## How to read support and coordinates

DNA evidence and RNA evidence answer different questions. Original DNA VCF
records are retained for nine SV geometries; 14 have RNA-only nominations.
Small-variant DNA depth/VAF comes from the frozen catalogue, with each library
kept separately. It has not been recounted here. A missing DNA call or unavailable
VAF means **unknown**, not DNA-negative. Processing copies are not added.

For ESVEE, **VF** is the caller's count of supporting variant fragments,
**SF** counts fragments with a read spanning the breakend, and **DF** counts
discordant fragments. These are breakend-level DNA counts, not RNA counts or
protein expression. PURPLE allele fraction/JCN and LINX cluster labels remain
caller estimates. An INV-labelled adjacency alone does not reconstruct an entire
derivative chromosome. Definitions follow the
[HMF VCF tags](https://github.com/hartwigmedical/hmftools/blob/master/hmf-common/src/main/java/com/hartwig/hmftools/common/sv/SvVcfTags.java).

RNA **alternate-allele fragments** can exceed **complete protein-window
fragments**. RNA fragments use read-group/query-name identities; independent
RNA molecule abundance is unresolved. Q20 requires recorded quality ≥20 at every base of that exact
window. SV quality was not assessed across the ORF interval; PacBio qualities
are unavailable. Barcode counts are labels, not verified malignant cells.
Counts are kept separate for bulk SARC0277 T2, tagged ONT T2, Kamil 10x T2
and PacBio T1. They are not summed across products, timepoints, orientations,
alternative starts or overlapping windows. A context-derived candidate with
zero complete witnesses is explicitly retained as a hypothesis.

Original VCF POS/REF/ALT are displayed unchanged. Exon and RNA-junction loci
in prose are **GRCh38 one-based**. Nominated SV cuts and insertion positions
are explicitly **interbase**. JSON genomic/cDNA intervals and mutation offsets
are **zero-based, half-open**. Removing an indel's VCF anchor is not HGVS
right normalization. These conventions follow the
[VCF specification](https://samtools.github.io/hts-specs/VCFv4.5.pdf).
RNA alignment junctions can differ from nominated DNA cuts because of splicing,
placement ambiguity or unplaced bases; the event sheets preserve both.

## SVs: original event and what happens to the RNA

| Event and complete sheet | Original DNA event/support | Transcript and protein interpretation |
| --- | --- | --- |
| [TPST1–CRCP / SV0575](data/sid-neoorf-event-report/events/adj-13d954785c89a21a1e8e937b.md) | **DNA origin unknown** in this selected entry. RNA nominations join the chr7:66205522 donor region back to chr7:66127704; an RNA caller's duplication label is not an original DNA call. | Retains TPST1 exon-1/5′-UTR sequence before the normal CDS and joins the downstream RNA to the CRCP-region locus. The candidate ATG is chr7:66205483 in ENST00000304842, 141 nt before its annotated CDS start. A 30-aa ORF changes after a shared 13-aa prefix relative to a potential reference uORF. **172 ONT T2 fragments / 139 labels; 15 PacBio T1 / 15 labels.** Initiation and RNA polarity remain unresolved. |
| [FOXO3–STRADA/CCDC47 / SV0241](data/sid-neoorf-event-report/events/adj-191e5810b18446510e1bbd6b.md) | Interchromosomal BND: chr6:108611245 `G` → `G[chr17:63745980[`. DNA VF: **T1 7, T1 organoid 5, T2 14**; original calls linked in sheet. | FOXO3-side cut is intronic; the chr17 locus overlaps STRADA intronic and CCDC47 UTR models on the other strand. RNA supports exploratory ORFs, including a **34-aa plus-stop** sequence with **ONT T2 5 fragments / 5 labels** and **PacBio T1 3 / 3**. Alternative 36- and 67-aa starts and other frames are retained. Reverse-complement-only support and unresolved polarity/frame prevent a canonical fusion-protein assignment. |
| [PARD3B–CDKN2B-AS1/CDKN2B, geometry A / SV0001](data/sid-neoorf-event-report/events/adj-5a298314b0272623da8eb76a.md) | BND: chr2:205310103 `G` → `GTGTGCTGATATG[chr9:22007648[`. DNA VF: **T0 18, T1 19, T2 15**. | PARD3B intronic region joins a chr9 region containing CDKN2B intronic and CDKN2B-AS1 noncoding models. Observed RNA placements include chr2:205310101(+) → chr9:22007648(+) with 13 unplaced bases `GTGTGCTGATATG`. The **10-aa plus-stop** hypothesis has **ONT 11 / 11 labels; PacBio 5 / 5**, with unresolved coding anchor/start. |
| [PARD3B geometry B / SV0319](data/sid-neoorf-event-report/events/adj-3834e653a87b6d5ec9b7b8a6.md) | T1-organoid BND at chr2:205310101 `T` → `TGTGTGCTGATATG[chr9:22007648[`; DNA **VF 6**, LINX complex-cluster annotation. | Slightly different nominated chr2 cut, with the same representative short RNA ORF and reused witnesses as geometry A. **Do not add their support.** Raw alleles are preserved rather than asserting identical DNA insertions. |
| [GABBR1–SLC29A1 / SV0111](data/sid-neoorf-event-report/events/adj-52ace853d3c63e505e9e5bba.md) | INV-labelled adjacency: chr6:29612906 `C` → `[chr6:44218908[C`; DNA VF **T1 12, T1 organoid 8, T2 10**. | GABBR1 endpoint is intronic; the nominated chr6:44218907 interbase partner has no overlapping Ensembl 115 transcript at the retained base. The **67-aa plus-stop PacBio hypothesis has 9 fragments / 8 labels**, but is only event-compatible. Its representative RNA join is chr6:44230858(−) → chr6:44230810(+) with `CCCGCCT`, not a direct bridge between the nominal DNA endpoints. Other ORFs include an ONT breakpoint-linked window; no resolved GABBR1–SLC29A1 protein is established. |
| [OTUD7A–FMN1 / SV0432](data/sid-neoorf-event-report/events/adj-79fa5cd408c4e72b464abb46.md) | T2 INV-labelled adjacency: chr15:31743670 `G` → `[chr15:33043479[CTG`; DNA **VF 9**. | Cuts overlap intronic OTUD7A and FMN1 models. A representative RNA join runs chr15:33043479(−) → chr15:31743670(+) with unplaced `CT`. The **21-aa plus-stop** hypothesis has **9 ONT fragments / 7 labels**. No resolved coding anchor or full fusion transcript is established. |
| [::IGKC dense control](data/sid-neoorf-event-report/events/adj-041dbfbb2239d183a5982838.md) | **DNA origin unknown**; RNA-only nomination at chr2 interbase cuts 88857683 and 88861220. | The **27-aa plus-stop 10x** hypothesis: **279 10x fragments / 12 labels**. Its RNA N-gap joins chr2:88857683(+) → chr2:88861886(+), differing from the nominated partner cut. Splice ambiguity, unresolved coding anchor/polarity and immune-receptor context prevent a tumor-specific neoORF claim. |
| [::CNN2 dense control](data/sid-neoorf-event-report/events/adj-02ce630f08dfbac9b6c4cde5.md) | **DNA origin unknown**; nominated chr19:1038882-right / chr9:42971353-right cuts. | **No resolved ORF** in the two acquired 10x views. The sheet gives retained-base transcript anatomy, original RNA nominations and bounded outcomes. Altered transcript structure, reading frame and protein sequence remain unknown; this is not a gene-expression negative. |
| [ATP5MG–KMT2A](data/sid-neoorf-event-report/events/adj-e7cac6d424ea83935654d789.md) | **DNA origin unknown**; supplemental Personalis T0 Isofox RNA nomination. | Observed splice geometry includes chr11:118401717(+) → chr11:118471662(+), compatible with ordinary read-through splicing. **29 context-derived hypotheses, zero complete-ORF witnesses.** Long KMT2A-containing sequences are provided in full, explicitly as unsupported complete-ORF hypotheses. No deletion or translated fusion protein is inferred. |

Representative sequences are below. The sheets and FASTA contain **all**
alternative sequences, placements, frames and support rows.

```text
TPST1–CRCP: MPSRRRGGSSLIPGRMPILRFSVTTRYFSY*
FOXO3 representative: MEKLCSLYYGCVPLLLRYISITEYPYNNGILLSN*
PARD3B representative: MNIRNVHFIS*
OTUD7A representative: MAAKAGTSLEARSLRPAWPTC*
GABBR1-region PacBio hypothesis:
MRSNLVSCPAWNSTATTSSSSLKDPGSRRPSWTSLAKERSQEQAKRNLEFQSPTLSPPMKATLSKPS*
::IGKC control: MKTDGAATVRLISTLVPWPNVHHLLYC*
```

IGKC is the antibody kappa constant-region gene in the chromosome-2
[B-cell receptor locus](https://www.imgt.org/IMGTrepertoire/LocusGenes/locusdesc/human/IGK/Hu_IGKdesc.html).
Its 27-aa computational candidate is antisense to the supplied IGKC models
and has no annotated coding start. The 279 complete fragments / 12 barcode
labels come from T2 10x; the nomination and additional ONT input are also T2.
Neither longitudinal recurrence nor malignant-cell assignment is
established. This entry remains an immune-receptor control, rather than a
prioritized tumor neoORF. A normal BCR origin for this exact short antisense
ORF is also unproven. The
[anatomy ledger](data/sid-sv-recovery/fth-igkc-candidate-anatomy-v2.json)
retains the start assessments, support and source claims.

For TPST1, an example RNA join has eight unplaced bases `CAGGAAGG` and
alignment endpoints chr7:66205519(+) → chr7:66127709(+). These are retained
as RNA placements, not rewritten as a DNA insertion. Its potential reference
39-aa upstream ORF is `MPSRRRGGSSLIPDVGYLSEVDCPWPEHFPKIILSKISV`; the normal
TPST1 coding start is downstream and not retained as the start of this short
candidate. The [earlier annotation analysis](../figures/osteosarc/neo-orfs-and-long-reads.md)
documents that comparison. Regional Ribo-seq in
[Kowar et al.](https://doi.org/10.1038/s41467-025-63405-2) does not establish
initiation at this internal ATG or translation of Sid's altered sequence.

## Recovered original-pilot targets

The [targeted recovery](sid-sv-dense-recovery.md) adds the original pilot's
22 acquisition-failed source/geometry pairs, preserving both orientations and
the full original denominator. The expanded sheets describe all 14 additional
SV geometries, including entries without a resolved RNA path or protein.
Source/product counts remain separate, and partner-recovery limits remain
attached to the ONT results.

The recovered **FTH1–FTL** candidate joins minus-strand FTH1 on chromosome 11
to plus-strand FTL on chromosome 19. It starts at the annotated FTH1 ATG in
ENST00000273550, cDNA offset 209 (zero-based), with genomic start-codon bases
chr11:61967425, 61967424 and 61967423 in translation order (one-based).
Its entire 540-nt ORF matches **FTH1 codons 1–17 followed in frame by FTL
codons 14–175 and the normal FTL stop**, using ENST00000331825. The resulting
179-aa sequence is:

```text
MTTASTSQVRQNYHQDSEAAVNSLVNLYLQASYTYLSLGFYFDRDDVALEGVSHFFRELAEEKREGYERLLKMQNQRGGRALFQDIKKPAEDEWGKTPDAMKAAMALEKKLNQALLDLHALGSARTDPHLCDFLETHFLDEEVKLIKKMGDHLTNLHRLGGPEAGLGEYLFERLTLKHD*
```

The first 20 residues also match FTH1 because `EAA` is shared with the FTL
segment. This describes sequence composition, not a uniquely placed genomic
fusion boundary. The representative RNA alignment has 63 unplaced bases
between chr11:61967374(−) and chr19:48965769(+) (zero-based boundaries).
They encode residues 18–38 and match FTL; a DNA insertion is not established.
The [exact comparison and annotation ledger](data/sid-sv-recovery/fth-igkc-candidate-anatomy-v2.json)
shows why the candidate has an annotated FTH1 start while its aggregate
start/frame labels remain ambiguous across overlapping transcript models.
**Two complete ONT T2 fragments** support this exact ORF; independent molecule
lineage, initiation and translation remain unresolved or unobserved. There is
no linked original DNA call and no established recurrence of this exact ORF
across timepoints. The historical 45-read T2 fusion-caller nomination measures
different evidence and is not a count of complete witnesses for this sequence.

The [published CTAT-LR-FUSION results](https://osteosarc.com/ctat_lr_fusion/FTH1--FTL/)
already call this gene pair in ONT T2, with 45 reads at the main nominated
junction and two other FTH1–FTL junctions supported by 18 and 12 reads. The
main call uses canonical splice boundaries and carries TCGA/CCLE catalogue
tags. In the matching published T2 concordance table, Isofox has no gene-pair
match, LINX has no DNA-driven gene-pair fusion, and PURPLE has no DNA SV
breakend within 10 kb of either nominated RNA breakpoint. These are the
reported caller/filter/window results, not a proof of DNA absence.
The two original query names supporting this exact 179-aa ORF are absent from
the exported CTAT supporting-read lists; processing lineage remains unresolved.
Thus CTAT corroborates the gene pair, but does not establish 45 complete
witnesses for this exact protein sequence.

As of October 8, 2026, this fusion is absent from the available public
[vaccine target records](https://osteosarc.com/vaccines/). A separate comparison
checked the current 183 source-variant entries, 37 overlap rows, and 95 recorded
vaccine-peptide sequences: no FTH1/FTL target or exact candidate peptide match
was found. None of the recorded peptides contains the candidate's 22
8–11-aa junction windows absent from both supplied normal parent proteins.
These parent-specific windows are not a whole-proteome novelty or presentation
assessment. The
[vaccine/caller comparison ledger](data/sid-sv-recovery/fth-vaccine-caller-comparison.json)
preserves the source hashes, exact matching rules, caller read IDs and results;
unprovided manufacturing construct sequences remain outside this check.

A [separate follow-up](sid-fth-fusion-followup.md) records the new junction
prediction screen and longitudinal comparison. For this exact 179-aa isoform,
NetMHCpan 4.2c found no junction-spanning 8–11-aa peptide below its default
class-I EL rank threshold for the five documented patient alleles. NetMHCIIpan
4.3e predicted five overlapping strong DRB3*01:01 15-mers sharing the core
`YHQDSEAAV` (best EL rank 0.07%); those overlapping windows are not independent
epitopes. These model-specific results do not establish presentation or explain
historical vaccine selection. No explicit rejection reason was found.
The preserved filtered CTAT table has 16 exact junctions shared across time
points (13 gene pairs); FTH1–FTL itself is called only at T2. T3 already has 24
filtered RNA calls, despite the fusion page's pending-mapping placeholder.


Four added events have original DNA calls but no resolved ORF in the recovered
10x views:

| Complete sheet | Original DNA event and available transcript anatomy |
| --- | --- |
| [SV0096](data/sid-neoorf-event-report/events/adj-01f832abce6755e57a037305.md) | Personalis BND: chr4:74716515 `G` → `G[chr11:54805939[`, VF 4. Neither retained endpoint overlaps a supplied Ensembl transcript. |
| [CNTNAP2 / SV0453](data/sid-neoorf-event-report/events/adj-02d0c8701be92d25ba193912.md) | T2 deletion adjacency: chr7:146504421 `T` → `T[chr7:146504535[`, VF 7. Both endpoints lie inside the same CNTNAP2 intron; the genomic deletion does not establish a translated frameshift. |
| [SV0271](data/sid-neoorf-event-report/events/adj-0323090002a789fe5398004e.md) | T1 deletion adjacency with unplaced `AT`: chr4:160600455 `A` → `AAT[chr4:160600600[`, VF 5. No overlapping supplied transcript at the endpoints. |
| [SV0406](data/sid-neoorf-event-report/events/adj-035949d7e8f1561d5e9bf526.md) | T2 BND: chr4:15847377 `G` → `[chr20:31245990[GG`, VF 13. The chr4 endpoint is intronic in CD38 and opposite-strand models; the chr20 retained base has no supplied transcript assignment. |

New exploratory 10x sequences include `MKQISRSLTLPCTS*` at
[TRDV1–TRAC](data/sid-neoorf-event-report/events/adj-020512ee5df1e1d9d00edcba.md)
(one complete fragment / one barcode label), `MQQNHPEVTAAMALTDIDLQLQF` at
[UBE2Q2P6–CSPG4P10](data/sid-neoorf-event-report/events/adj-02f752de4e9fd2a55e4924f5.md)
(one / one), `MILWPSNHPSKCLGKKR*` at
[CDK11A–CDK11B](data/sid-neoorf-event-report/events/adj-03b9fc802bfa90900b74a750.md)
(two / one), and `MQFLHKLILRENNLTELHKDSFEGL` at
[LRRC37A–LRRC37A2](data/sid-neoorf-event-report/events/adj-03db542e1425ff6d0005c1ab.md)
(three / two). Joins are splice ambiguous, frames are unresolved or ambiguous,
and the original DNA events are unknown. These observations do not establish
tumor-specific fusion proteins. The complete sheets give genomic RNA blocks,
reference transcript placement and every alternative sequence/support row.

## Frameshifts: DNA allele, local transcript change and sequence

| Complete event sheet | Original DNA allele and local RNA effect | Protein window and exact support |
| --- | --- | --- |
| [SPAG1](data/sid-neoorf-event-report/events/SPAG1-chr8-100213209.md) | chr8:100213209 `CG` → `C`: removes G at chr8:100213210 on the + strand. The event is exonic/CDS in supported models; ENST00000388798 places it in exon 11/19. Original DNA metadata includes T1 UCLA depth 96/VAF 0.1354, T1 organoid 72/0.125 and T2 UCLA 111/0.027. | 10x T2: **20 alternate-allele fragments**. `GAPQRARPRRPART`: **16 complete / 15 Q20**, 10 shifted aa. Longer `GAPQRARPRRPARTSGAHGGPLRRRRRAA`: **2 complete / 2 Q20**, 25 shifted aa beginning at offset 4. No stop is observed in these windows; full altered protein remains unknown. |
| [TECPR1](data/sid-neoorf-event-report/events/TECPR1-chr7-98241127.md) | chr7:98241127 `TG` → `T`: deletes genomic G at chr7:98241128, corresponding to C in the −-strand transcript. Source NM015395 label is `c.774del` / `p.Thr259fs`; transcript-specific exon/CDS placement is in the sheet. DNA T0 depth 24/VAF 0.25; T1 UCLA 74/0.0946; T2 UCLA 157/0.0446. | Bulk T2: **17 alternate-allele fragments**. `HDLLWAHSGRDRPWSG`: **15 complete / 14 Q20**, 11 shifted aa. Longer `PHDLLWAHSGRDRPWSGKESTGAIPKEVPGPSWSLLDLKTGSCTS`: **2 complete / 2 Q20**, 39 shifted aa beginning at offset 6. No in-window stop; full altered protein unknown. |
| [PIP5K1A](data/sid-neoorf-event-report/events/PIP5K1A-chr1-151242178.md) | chr1:151242178 `AG` → `A`: deletes G at chr1:151242179 in + strand coding models. Source `p.Gly474fs`; Ensembl transcript anatomy is retained separately. DNA T1 UCLA depth 104/VAF 0.1635; T1 organoid 78/0.1538; T2 UCLA 165/0.103. | ONT T2: **7 alternate-allele fragments; 2 complete / 1 Q20** cover the 81-aa window below, including **21 shifted aa plus stop**. 10x T2: **2 complete / 2 Q20** cover `SRRAAPVATPALLTSHRSLGN`, with 17 shifted aa and no observed stop. Bulk has shorter windows. These counts are not pooled. |
| [KTN1](data/sid-neoorf-event-report/events/KTN1-chr14-55627965.md) | chr14:55627965 `G` → `GTT`: inserts TT after the anchor, at interbase 55627965. Exonic/CDS in compatible + strand models. Original DNA depth/VAF unavailable. | ONT T2: **11 complete / 8 Q20** for `TTLIHQLQEKDKFYSLL*`; the last five aa `FYSLL` are shifted. A short altered tail with a linked stop, not evidence for the full-length protein or its initiation. |

The full PIP5K1A ONT coding window is:

```text
KLEHSWKALVHDGDTVSVHRPGFYAERFQRFMCNTVFKKIPLKPSPSKKFRSGSSFSRRAAPVATPALLTSHRSLGNTRHK*
```

SPAG1's source explicitly warns that its `p.Arg406fs` label differs from
the linked gnomAD transcript annotations (`p.Gly407AlafsTer42`). The report
preserves both rather than assigning a universal amino-acid position. The
observed shifted sequence, exact window mutation offsets and Ensembl transcript
IDs are authoritative for the windows reported here; source RefSeq labels,
genomic padding and HGVS normalization can have different placements.

## Unresolved events and coding controls

All remaining nominations have event sheets, even where no protein is resolved:
**ACSL6, ANKRD20A8P, BAZ1B, CABLES1, CCDC40, CDHR18P, CDKN2B-AS1,
CGNL1, COL3A1-Splice, DCHS2, FAM157A, FGFR3, GAPVD1, both GLIS3 alleles,
GOLGA6L2, both MAP2 alleles, MUC3A, MYO15B, the unassigned chr15:77027489
event, NNT, PCDH19, SULT6B1 and TAF7L**. Each sheet describes the event
location in every intersecting gene's transcript models, including intronic,
UTR, CDS, mixed exon/intron and outside-isoform placements, original alleles
and DNA/RNA support. The [index](data/sid-neoorf-event-report/README.md)
links each separately.

Some have alternate-allele RNA without a resolved protein. For example,
ANKRD20A8P has 19 bulk / 99 ONT / 17 10x / 8 PacBio alternate fragments;
PCDH19 has 5 bulk / 6 ONT / 1 10x / 0 PacBio. Those counts do not establish
an altered mature splice transcript or a neoORF. Splice-site substitutions
were screened as focal alleles; arbitrary splice changes, retained introns
and exonization were not comprehensively discovered. COL3A1's long genomic
deletion crosses intronic/exonic annotation; its genomic length cannot be
read directly as a translated frameshift.

MUC3A ALT is the string `dup`, with source duplication/cDNA labels but no
literal sequence. It is **unassessable in all four products**, not zero RNA
support. Its sheet describes the nominated locus; duplication boundaries,
altered transcript and protein remain unresolved.

**RNF213, H1_2 and GTF3C5** have in-frame coding windows. Their original
DNA deletions, RNA support and all 774 source-specific windows are included,
explicitly marked as controls rather than frameshift neoORFs.

## Limits, verification and reproduction

The selected search has 100 oriented SV views and 124 completed literal
allele/product screens, plus four unassessable MUC3A outcomes. Every frozen
sequence/support row is retained. The generator checks input hashes and verifies
each SV ORF nucleotide interval against its original RNA reconstruction before
mapping the interval to genomic blocks, on both strands. Existing
[native coding-window checks](data/sid-neoorf/coding-window-checks.json.gz)
independently verify nine selected frameshift windows against original SAM
sequence, fragment identities and recorded qualities.

Protein novelty is measured only against the pinned Ensembl 115 annotated
proteome (245,535 entries; exact and I/L-collapsed 8-mers). It does not establish
absence from normal noncoding ORFs, matched normal tissue, or other patient
alleles. RNA sequence support does not establish initiation, translation,
tumor specificity, HLA binding or presentation. Boundaries and truncation flags
remain in the ledger; the full 1,493-geometry catalogue is outside this focused
report, not a negative result.

With the frozen selected input directories and pinned Ensembl 115 installed:

In current releases, `isovar[data]` installs the Biopython verification
dependency; it contains no Sid data and does not download samples. Original
read acquisition is owned by Osteosarc through the
[selected-read workflow](sid-neoorf-read-subsets.md); the report FASTA, support
table and downloadable ledger are linked above. Simplifying the install and
migrating the checks to Biotite is tracked in
[#579](https://github.com/openvax/isovar/issues/579). The
original-BAM/SAM witness checkers also require the **SAMtools CLI** on `PATH`;
CI installs it explicitly. The report generator itself reads frozen JSON and
annotation, without invoking SAMtools or downloading whole libraries.

```sh
python -m examples.report_sid_neoorf_events OUTPUT \
  --screen notes/data/sid-neoorf/candidate-screen.recovered-final.json.gz \
  --checks notes/data/sid-neoorf/coding-window-checks.json.gz \
  --audit PRIOR_AUDIT --audit ADDITIONS_AUDIT --audit RECOVERY_AUDIT
```

Use a fresh output directory; audit JSON is immutable. The
[selected read subsets](sid-neoorf-read-subsets.md) and
[screen methods](sid-sv-neoorf-candidates.md) describe input scope and acquisition
limits. The event report adds annotation and provenance to those completed
screens; it performs no new whole-library download.
