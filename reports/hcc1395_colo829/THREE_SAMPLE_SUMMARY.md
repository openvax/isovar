# Sid, COLO829 and HCC1395: variant funnels and supported ORFs

Evidence reviewed 2026-10-26. This index joins the existing cell-line atlases,
Sid's selected-event report and the new shared COLO829 SV workflow. It does not
claim that every sample was rerun with the same software version. RNA sequence
support does not establish translation, tumor specificity or antigen presentation.

## Reports and complete results

| Sample | Report | Per-event results |
|---|---|---|
| Sid | [NeoORF report](../../notes/sid-neoorf-event-report.md) | [55 event sheets](../../notes/data/sid-neoorf-event-report/README.md), [all support rows](../../notes/data/sid-neoorf-event-report/support.tsv), [protein/window FASTA](../../notes/data/sid-neoorf-event-report/proteins.fasta) |
| COLO829 | [Whole-run atlas](COLO829_ATLAS.md) | [46,064 small-variant rows](colo829_atlas_calls.tsv.gz), [local proteins](colo829_atlas_proteins.tsv.gz), [new SV inventory](unified_sv_audit/REPORT.md), [SV3](unified_sv_audit/SV3.md) |
| HCC1395 | [Whole-run atlas](HCC1395_ATLAS.md) | [41,072 small-variant rows](hcc1395_direct_rna_atlas_calls.tsv.gz), [local proteins](hcc1395_direct_rna_atlas_proteins.tsv.gz), [SV screen](hcc1395_direct_rna_sv_screen.tsv), [EIF3K–CYP39A1 fusion](eif3k_cyp39a1/report.md) |

The cell-line funnels account for every row in their nominated DNA inventories.
Sid's detailed report covers **23 selected SV geometries and 32 small-variant
nominations**, not all Sid variants. The broader frozen Sid SV catalogue has
2,258 nominations representing 1,493 geometries. Its [four-platform pilot](../../notes/sid-sv-four-platform-pilot.md)
covers 50 geometries/400 oriented views, followed by selected recovery. A full
catalogue-wide result remains missing. Failed/unacquired inputs remain unknown.

## Small-variant funnels

These source-record counts use the cell-line atlases' pinned Ensembl 87 annotation.
The exon and all-locus funnels are separate views of the same inventory.

| Stage | HCC1395 | COLO829 |
|---|---:|---:|
| All small-variant rows | 41,072 | 46,064 |
| SNVs / indels | 39,447 / 1,625 | 44,005 / 2,059 |
| Any annotated exon involved | 1,728 | 1,395 |
| Exonic and RNA ALT observed | 292 | 502 |
| Exonic and local protein context | 22 | 124 |
| RNA ALT observed, all loci | 343 | 1,933 |
| ALT-observed calls with local protein context | 22 | 124 |
| Ranked local protein contexts | 119 | 12,218 |
| DNA-predicted top-effect frameshifts | 9 | 5 |
| ncORF coding-block overlap | 17 | 5 |
| ncORF-overlap calls with RNA ALT | 5 | 2 |
| Reconstructed mutant ncORF proteins | 0 | 0 |

Local-context totals can include reference-identical/synonymous products; they
are not counts of confirmed mutant proteins. Cohort-wide separation remains
tracked in [#551](https://github.com/openvax/isovar/issues/551).

| RNA status, all source rows | HCC1395 | COLO829 |
|---|---:|---:|
| Alternate observed | 343 | 1,933 |
| Reference only | 239 | 1,637 |
| Reference and other only | 53 | 36 |
| Other allele only | 19 | 27 |
| Alignment without allele call | 7,300 | 166 |
| No alignment overlap | 33,118 | 42,265 |

HCC1395 uses ONT direct RNA SRR23456709 (1,049,452 reads); COLO829 uses Illumina
SRR8615617 (115,269,780 pairs), aligned in full before nominated-site selection.
COLO829's genomic `minimap2 -ax sr` alignment can miss splice-crossing reads.
RNA and benchmark DNA cultures differ. Missing support or an RNA reference allele
does not establish absence of the alteration in its source DNA culture.

## Strongest frameshift windows

ALT support differs from support for an entire coding window. None of these
local windows resolves a full-length mutant gene protein.

| Sample / event | ALT fragments | Complete-window fragments | Resolved sequence |
|---|---:|---:|---|
| COLO829 CDKN2A, chr9:21971154 CCG>C | 6,967 | 5,248 (5,338 allele-read observations) | Rank-1 `MMMGSARVAELLLLHGGAQLRRPR`; mutation interval [16,24). DNA top effect p.A68fs; compatible transcript isoforms retained. |
| COLO829 FHL2, chr2:105361307 G>GCCCC | 845 | 654 (747 allele-read observations) | Rank-1 `FNCKKCSLSLVGRGFLTERDDILWGPRLWERHLNS`; mutation interval [23,35). DNA top effect p.C388fs. |
| Sid SPAG1 | 20, 10x T2 | 16; 15 Q20 | 10 shifted aa; longer 25-aa shifted tail has 2 complete/Q20 fragments. No stop observed. |
| Sid TECPR1 | 17, bulk T2 | 15; 14 Q20 | 11 shifted aa; longer 39-aa shifted tail has 2 complete/Q20 fragments. No stop observed. |
| Sid PIP5K1A | 7, ONT T2 | 2; 1 Q20 | 81-aa local window includes 21 shifted aa plus stop. Separately, 10x has 2 Q20 witnesses for a shorter window. |
| Sid KTN1 | See event sheet | 11; 8 Q20, ONT T2 | Five shifted aa plus stop. |
| HCC1395 C8orf76, chr8:123220221 A>AGGGCTAC | 3 | 2 | Local frameshift context; not a high-coverage lead. |

Other COLO829 predicted frameshifts have 3, 1 and 0 ALT fragments; only the
three-fragment call has local contexts. HCC1395's other eight have no ALT
observation. Support for a short window cannot be assigned to a longer tail.

## ncORFs and SV-derived candidates

| Candidate | RNA support, kept separate by source | Interpretation |
|---|---|---|
| Sid TPST1–CRCP altered uORF | 172 complete ONT T2 fragments; separately 15 PacBio T1 | 30 aa plus stop: 13-aa shared prefix, 17-aa altered tail. Candidate ATG is 141 nt upstream of TPST1's annotated CDS. Initiation, polarity and exact somatic DNA origin unresolved. |
| Sid FOXO3–STRADA/CCDC47 | Representative 34-aa-plus-stop ORF: 5 ONT and 3 PacBio fragments | Original DNA BND exists; alternative starts/frames and reverse-complement-only RNA interpretation remain unresolved. |
| Sid PARD3B–CDKN2B-AS1/CDKN2B | 10-aa-plus-stop: 11 ONT and 5 PacBio | No established coding start. Two nominated geometries reuse witnesses. |
| Sid OTUD7A–FMN1 | 21-aa-plus-stop: 9 ONT fragments | No established coding anchor. |
| Sid FTH1–FTL | 2 complete ONT T2 fragments | 179-aa in-frame candidate with annotated FTH1 start and FTL stop. Exact DNA origin/recurrence unestablished. |
| HCC1395 EIF3K–CYP39A1 | 149 Illumina junction fragments; exact coding tiling by 19 mapped reads; separately 9 ONT coding-span reads, 4 with exact junction | 531-aa complete annotated-frame candidate. Illumina input is fusion-enriched; counts are not whole-sample frequency. Source provenance is in the dedicated report. |
| COLO829 SV3 | 3 directly placed junction fragments; 5 exact 43-base junction-window fragments | Exon-linked sequence 280 nt downstream does not phase the join to CDS. The 37/162-aa exploratory ORFs remain unframed; no justified annotated protein. |

Original HCC1395/COLO829 benchmark SV funnels, respectively: **1,788/188 records
→ 354/97 exon-involved → 172/64 with RNA near an end → 6/6 with split/mate
candidates → 3/0 with junctions in selected exonic replay → 0/0 resolved proteins**.
This is a bounded benchmark screen; the independently supplied HCC1395 fusion
is a separate positive control. COLO829 SV3 is outside that exon-involved rule;
its low-MAPQ joins were missed by the MAPQ-20 screen. The new [inventory report](unified_sv_audit/REPORT.md)
provides all retained-path outcomes and explicit unresolved inputs.

The new COLO829 run covers **188 source records / 380 source-orientation runs**:
266 no-candidate, 66 regional-only, 30 event-linked-candidate and 18 unassessed
outcomes. All **657 retained paths** are assessed. Eight paths across eight
runs have direct-junction support, attributable to three source records:
SV3 and the older `truthset_5_1`/`truthset_5_2` breakends. Their nominated DNA
coordinates differ, but the recovered RNA junction witnesses are reused; these
are not independent RNA replications. There are **zero event-linked
annotated-frame protein candidates**. A reconstruction's event-linked label
does not alone establish its frame or distinguish it from competing explanations.

ncORF overlap is annotation, not translation evidence. The strongest overlapping
calls are COLO829 chr17:58089866 G>A (dORF, 45 ALT fragments) and HCC1395
chr9:137255474 G>C (intORF, 22). The HCC call also has a canonical coding context,
which does not establish translation of its overlapping ncORF.

## Threshold: at least 20 supporting fragments

Recounted from the linked call/protein/support tables, without pooling libraries,
orientations or overlapping alternatives. Fragments are template identities,
not UMI-confirmed independent molecules.

| Scope | At least 20 ALT fragments | At least 20 complete reported-window fragments |
|---|---:|---:|
| COLO829, all small variants | 220 calls | 74 calls, including non-frameshift products |
| HCC1395, all small variants | 30 calls | 0 calls |
| COLO829, predicted frameshifts | 2: CDKN2A and FHL2 | The same 2 calls |
| HCC1395, predicted frameshifts | 0 | 0 |
| Sid, selected frameshift windows | SPAG1 has 20 ALT fragments | 0 windows |
| Sid, selected SV ORFs | Not a focal ALT count | 2 event groups: TPST1–CRCP and IGKC |

Sid's threshold gives 86 source-specific SV-ORF rows / 85 event-specific sequence
groups, all belonging to **two event groups**. There is **one additional group
beyond TPST1–CRCP**, not 85 independent high-support tumor ORFs. IGKC is an
immune-receptor, splice/polarity-ambiguous control; its 27-aa 10x candidate has
279 complete fragments. Its alternatives are not separate tumor discoveries.
Two in-frame control events, H1-2 and GTF3C5, also pass but are not frameshift/SV
neoORFs. No library counts or Q20 and quality-unassessed counts are pooled.

## Remaining coverage

The shared workflow separates junction, branch, exon/CDS, frame and whole-sequence
support. Original COLO829 PacBio tumor/normal RNA remains unavailable ([#574](https://github.com/openvax/isovar/issues/574));
unresolved SV definitions remain explicit ([#575](https://github.com/openvax/isovar/issues/575));
SV3 biological frame resolution stays open ([#563](https://github.com/openvax/isovar/issues/563)).
The new workflow has not yet covered the complete Sid and HCC1395 inventories.
Sid uses Ensembl 115 and the cell-line atlases use Ensembl 87, so reference-novel
peptide totals are not comparable. Exact matching undercounts noisy reads;
missing data and bounded searches are not evidence of absence.
