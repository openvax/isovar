# Remaining Sid pilot SV read recovery

[PR #460](https://github.com/openvax/isovar/pull/460) accounts for the original
pilot's 22 acquisition-failed source/geometry pairs (44 orientations), using
osteosarc 0.15.5's inventory source-date checks. The original 400-view pilot
and its failures remain preserved. This is a selected recovery, not a full
1,493-geometry catalogue search.

All 22 selected read subsets are acquired. Of the 44 views, 42 completed
reconstruction; ONT IGKC forward timed out after 14,400 seconds and reverse
was interrupted by the coordinator. Their candidate/support counts remain
unknown. [Follow-up #470](https://github.com/openvax/isovar/issues/470) tracks
those two views. The [coverage ledger](data/sid-sv-recovery/recovery-coverage.json.gz),
[independent checks](data/sid-sv-recovery/recovery-checks.json.gz), and
[revised pilot](data/sid-sv-recovery/revised-pilot.json.gz) retain all 44 outcomes
and the other 356 historical pilot views. The full catalogue denominator
remains separate from these 400 observed pilot views.

## Completed 10x outcomes

Acquired records include requested regional context and any matched partners.
These are alignment records, not junction fragments or protein abundance.
Each completed view scans every preserved input record for candidate support;
discovery remains bounded at 50,000 selected records.

The public `isovar sv-rna --max-records` option can raise the raw discovery
budget in a separate attempt. Current selection gives priority to event-related
joins and groups by exact sequence, placement and quality; that strict grouping
does little to compress noisy long reads. In a
[read-only FTH sample](data/sid-sv-recovery/discovery-grouping-measurements.json) of 24,614
segment bundles, sequence/placement grouping left 23,665 classes, splice-pattern
grouping left 343, and splice plus local allele/coverage grouping left 3,685.
The last measurement took nine seconds and covered only 100 bases around the
breakends. These are exploratory grouping counts, not validated discovery
representatives or reconstruction timings. [#464](https://github.com/openvax/isovar/issues/464)
tracks preserving complete sequence/quality/lineage evidence while budgeting
distinct junction and allele patterns and benchmarking higher limits.

A separate [storage codec experiment](data/sid-sv-recovery/lossless-compression-measurements.json)
used deterministic 0.5% template samples, preserving every decoded alignment
field, quality and optional tag. Compared with level-6 BAM, self-contained
CRAM 3.1 archive output was 33.4% smaller for FTH (7,860 records) and 37.9%
smaller for IGKC (13,817 records); level-9 BAM saved only about 3.6%. These
are sampled storage ratios, not measured full-file savings or faster
reconstruction. The original BAM/index bytes remain unchanged. CRAM preserves
alignment content here but does not restore the original BAM byte checksum;
the current exact-byte audit/replay still uses its original frozen BAMs.

| Target | Acquired records | Forward / reverse outcome | ORF hypotheses |
| --- | ---: | --- | ---: |
| TMSB10::TMSB4X / shared:SV0581 | 4,204,441 | no_candidate_paths / no_candidate_paths | 0 |
| shared:SV0096 | 43 | no_candidate_paths / no_candidate_paths | 0 |
| TRDV1::TRAC / shared:SV0631 | 23,795 | splice_ambiguous_candidates / splice_ambiguous_candidates | 13 |
| ::CNN2 | 98,966 | no_candidate_paths / no_candidate_paths | 0 |
| CNTNAP2::CNTNAP2 / shared:SV0453 | 85 | no_candidate_paths / no_candidate_paths | 0 |
| EEF1AKMT2::FAM53B | 13,804 | splice_ambiguous_candidates / splice_ambiguous_candidates | 0 |
| UBE2Q2P6::CSPG4P10 | 858 | splice_ambiguous_candidates / splice_ambiguous_candidates | 6 |
| shared:SV0271 | 1 | no_candidate_paths / no_candidate_paths | 0 |
| shared:SV0406 | 26,403 | no_candidate_paths / no_candidate_paths | 0 |
| POLR2J2::POLR2J3 | 24,901 | no_candidate_paths / no_candidate_paths | 0 |
| ZFAND2A::C7orf50 | 25,184 | no_candidate_paths / no_candidate_paths | 0 |
| CDK11A::CDK11B | 35,301 | splice_ambiguous_candidates / splice_ambiguous_candidates | 70 |
| LRRC37A::LRRC37A2 | 21,983 | splice_ambiguous_candidates / splice_ambiguous_candidates | 16 |
| ::IGKC | 773,288 | event_linked_candidates / event_linked_candidates | 176 |
| FTH1::FTL / shared:SV0580 | 5,669,578 | no_candidate_paths / no_candidate_paths | 0 |
| UBC::RPS27A / shared:SV0616 | 1,351,544 | no_candidate_paths / no_candidate_paths | 0 |

TRDV1–TRAC, UBE2Q2P6–CSPG4P10, CDK11A–CDK11B and LRRC37A–LRRC37A2
produce 105 exploratory hypotheses across these 10x views. Their representative
full-interval counts are respectively 1, 1, 2 and 3 fragments. All have
splice-ambiguous joins and unresolved or ambiguous frames; the selected
nominations have no linked original DNA VCF call. The expanded event sheets
retain exact sequences, support, RNA blocks and reference transcript anatomy.

The ::IGKC control reproduces its previously recorded 27-aa plus-stop sequence,
`MKTDGAATVRLISTLVPWPNVHHLLYC`, with 279 complete fragments and 12 barcode labels.
Its junction is splice ambiguous and DNA origin/initiation remain unknown.
IGKC encodes the antibody kappa constant region in the chromosome-2 B-cell
receptor locus; normal V–J rearrangement and transcription with IGKC are
described by [IMGT](https://www.imgt.org/IMGTrepertoire/LocusGenes/locusdesc/human/IGK/Hu_IGKdesc.html).
This particular computational ORF runs on the plus strand, antisense to the
supplied minus-strand IGKC/IGKJ models, and has no annotated coding start.
Its counts do not establish a normal antibody protein or a tumor neoORF.
The nomination and completed 10x evidence are T2; the additional ONT input is
also T2. Longitudinal recurrence, malignant-cell assignment and tumor
specificity are unestablished. It remains an immune-receptor control and is
not prioritized as a tumor neoORF without that additional evidence. The
[candidate-anatomy ledger](data/sid-sv-recovery/fth-igkc-candidate-anatomy-v2.json)
preserves the supporting results and source claims.
No candidate path is resolved in the 10x TMSB10–TMSB4X, FTH1–FTL or UBC–RPS27A
views. A bounded search's zero-candidate result is not a gene-expression negative.

## ONT outcomes

| Target | Acquired records | Forward / reverse outcome | ORF hypotheses |
| --- | ---: | --- | ---: |
| TMSB10::TMSB4X / shared:SV0581 | 1,758,048 | event_linked_candidates / event_linked_candidates | 1 |
| ZFAND2A::C7orf50 | 5,640 | splice_ambiguous_candidates / splice_ambiguous_candidates | 1 |
| CDK11A::CDK11B | 12,747 | splice_ambiguous_candidates / splice_ambiguous_candidates | 100 |
| LRRC37A::LRRC37A2 | 7,635 | no_candidate_paths / no_candidate_paths | 0 |
| FTH1::FTL / shared:SV0580 | 1,567,470 | event_linked_candidates / event_linked_candidates | 37 |
| ::IGKC | 2,715,459 | reconstruction_timeout / reconstruction_interrupted | unknown |

These inputs retain incomplete partner-recovery flags. Each completed view
counts candidate support across every preserved alignment record, while
discovery in the dense TMSB input remains bounded at 50,000 records.

The ZFAND2A–C7orf50 reverse view has a 103-aa plus-stop hypothesis supported
by seven complete fragments and seven barcode labels:

```text
MEFPDLGKHCSEKTCKQLDFLPVKCDACKQDFCKDHFPYAAHKCPFAFQKDVHVPVCPLCNTPIPVKKGQIPDVVVGDHIDRDCDSHPGKKKEKDGKTEEKSS*
```

The TMSB10–TMSB4X reverse view has one complete fragment and one barcode label
for `MYSAYRITRLITIRLLASPVQSFLPTVGRVLLACVCLFFSLALSETVPCSFLYKFLKQSVQSRSEPEILVLPLALLQQE`,
without an observed stop. The CDK11 views produce 100 alternative hypotheses,
with at most two complete fragments for a sequence. Frames, initiation and
read lineage remain unresolved; barcode labels do not establish independent
cells or molecules. RNA nominations alone do not establish an original DNA event.

The FTH1–FTL forward view has 11 ORF hypotheses. Its strongest 179-aa
plus-stop sequence has two complete fragments and two barcode labels. The
reconstructed join connects the minus-strand chr11 boundary at 61,967,374
to the plus-strand chr19 boundary at 48,965,769, with 63 unplaced RNA bases
between them (zero-based genomic boundaries). Its ATG maps exactly to the
annotated FTH1 coding start in ENST00000273550 and ten other coding models:
chr11 bases 61,967,425, 61,967,424 and 61,967,423 in minus-strand translation
order (one-based), cDNA offset 209 in ENST00000273550 (zero-based).
The entire 540-nt ORF exactly equals the first 51 coding bases of that FTH1
transcript followed by ENST00000331825's FTL coding sequence from codon 14
through its normal `TAA` stop. Its protein is therefore sequence-compatible
with **FTH1 residues 1–17 → FTL residues 14–175 → stop**, in frame.
The aggregate start/frame labels remain ambiguous because overlapping
noncoding, intronic and antisense transcript interpretations are retained;
they do not mean that the FTH1 coding-start match is absent. Initiation and
translation remain unobserved. The
[candidate-anatomy checks](data/sid-sv-recovery/fth-igkc-candidate-anatomy-v2.json)
record both occurrences and the exact nucleotide/protein comparisons.
The nomination has no linked original DNA VCF call;
the event report retains its transcript models, exact RNA/protein sequences
and both sequence witnesses without assigning a DNA mechanism.
Its first 20 amino acids match the primary [FTH1 sequence](https://www.uniprot.org/uniprotkb/P02794/entry)
and its last 162 match [FTL](https://www.uniprot.org/uniprotkb/P02792/entry),
with three overlapping shared residues. This supplementary
[reference comparison](data/sid-sv-recovery/fth-protein-reference-comparison.json)
describes sequence composition; novelty screening uses the separately frozen
Ensembl 115 proteome.
The 63 unplaced junction bases encode residues 18–38 and match FTL sequence;
that sequence match does not resolve their genomic placement or establish a
DNA insertion. The two complete witnesses for this exact ORF are ONT T2.
Historical nomination metadata reported 45 T2 ONT fusion-caller reads, which
is a different support definition and nominated sequence. Earlier T0/T1
adjacency-screen entries report no observed join, while several T3 entries
remain incomplete; these historical entries are not exhaustive negative
screens. Recurrence of this exact ORF across timepoints is unestablished.

This is a conventional **in-frame fusion ORF lead**, with an annotated FTH1
start and normal FTL stop. [CTAT-LR-FUSION](https://osteosarc.com/ctat_lr_fusion/FTH1--FTL/)
already calls the FTH1–FTL gene pair in T2: 45 reads at the main nominated
splice-boundary junction, with two other junctions supported by 18 and 12.
Those caller counts are separate from the two exact complete-ORF witnesses,
whose original query names do not occur in the exported CTAT read lists.
The published T2 cross-check has no Isofox gene-pair match, no LINX fusion,
and no PURPLE DNA SV breakend within 10 kb of either main RNA breakpoint.
In an October 8 comparison of public vaccine records, none of the 95 recorded
peptide sequences matched the candidate or its 8–11-aa junction windows
absent from both normal parent proteins. No FTH1/FTL vaccine target was listed.
The [comparison ledger](data/sid-sv-recovery/fth-vaccine-caller-comparison.json)
preserves the live source snapshots and separates gene-pair caller evidence,
exact-ORF witnesses, and vaccine membership. It does not assess antigen
presentation or unprovided manufacturing constructs.

The [junction-prediction and longitudinal follow-up](sid-fth-fusion-followup.md)
evaluates this exact 179-aa isoform. NetMHCpan 4.2c found no junction-spanning
8–11-aa peptide below its default class-I EL threshold for the five documented
patient alleles. NetMHCIIpan 4.3e predicted five overlapping strong DRB3*01:01
15-mers with the same core, `YHQDSEAAV` (best EL rank 0.07%). This is a prediction
of one overlapping region, not five independent epitopes or observed
presentation. No explicit historical vaccine-rejection reason was found.
The filtered CTAT table has 16 exact junctions shared across time points,
representing 13 gene pairs; FTH1–FTL itself is called only at T2. T3 already has
24 filtered RNA calls despite the fusion page's pending-mapping placeholder.

## Ordinary Python/CLI reproduction

The example runner uses osteosarc's public indexed acquisition API and Isovar's
RNA reconstruction engine. It needs a checkout and frozen base inventory and
transcript-reference files; installed portable replay uses the public
[audit-bundle API and CLI](../docs/audit-bundles.md).

```sh
cd ISOVAR_CHECKOUT
python -m examples.recover_sid_sv_pilot RECOVERY \
  --base-audit FROZEN_BASE --pilot-coverage notes/data/sid-sv-pilot/coverage.json.gz \
  --cache CACHE --stage acquire --max-acquisition-records 6500000 \
  --seed-timeout 1200 --max-seconds 3600
python -m examples.recover_sid_sv_pilot RECOVERY \
  --base-audit FROZEN_BASE --pilot-coverage notes/data/sid-sv-pilot/coverage.json.gz \
  --cache CACHE --stage reconstruct --max-acquisition-records 6500000 \
  --seed-timeout 1200 --max-seconds 3600
python -m examples.check_sid_sv_recovery RECOVERY CHECKS
```

The checker requires every view to be reconstructed by default. To verify a
publication that explicitly retains recorded failures, as this one does, add
`--allow-unresolved`. Every selected input and outcome is still required;
the ledger reports completed and unresolved views separately.

The default free-space guard is 8 GiB. Partner recovery allows four rounds,
20,000 intervals, 20 million bases and 128 queries of at most 30 seconds each;
seed acquisition has its separate 1,200-second budget. Each reconstruction
view has a 3,600-second budget. Limits and failed queries remain in the receipts.
The declared ::IGKC retry uses a 3,600-second seed budget and a 7,200-second
budget per reconstruction view. Its seed contains 2,715,423 records, matching
an independent remote SAMtools regional count. The retry's streaming reader
uses SAMtools/HTSlib 1.24 through pysam 0.24.1; previous failed/interrupted
attempts remain preserved. These overrides are recorded separately before
the corresponding work starts.
Its final acquired subset includes 36 recovered partners in addition to that
complete regional seed; the partner-query limit remains explicit.
Changed scientific inputs/budgets require a separate attempt. The public runner
rejects a package imported from a different checkout before any preparation or
acquisition; use the same checkout for the example and Isovar package.

Coverage snapshots have content-derived names. Completed zero-candidate scans
are distinct from absent, unassessable or timed-out results, whose counts remain
unknown. Corrected reporting attempts reuse only checksum-verified receipts and
results with identical scientific engine, inputs and parameters. Their original
producing code and controls remain pinned; coverage history is retained.

Independent validation recounts total/regional records with native samtools,
translates every exported ORF with Biopython's standard genetic code, and checks
representative sequence/fragment witnesses against original SAM and BAM records.
[SAM orientation/clipping](https://samtools.github.io/hts-specs/SAMv1.pdf) and
[NCBI genetic code 1](https://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi#SG1)
define those checks. RNA hypotheses do not establish initiation, translation,
malignant-cell identity, tumor specificity or presentation.

The final independent ledger recounts all 22 input subsets and verifies all
44 recorded outcomes, including 42 completed reconstructions and 420 translated
ORF hypotheses. Its explicit unresolved mode does not report the two unfinished
IGKC views as verified reconstructions. Earlier failed budgets, the storage stop
and both original IGKC stop records remain pinned. FTH completed under a
separately declared 14,400-second per-view budget with unchanged scientific
settings. The original 3,600-second and proposed IGKC 7,200-second budgets
remain provenance.
