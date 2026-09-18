# Expanded RNA evidence and complete protein alternatives

Audited September 17, 2026. This extends the nine candidates from
`RNA_FOOTPRINTS.md` with eight additional indexed GRCh38 RNA products:
T0 BostonGene and Personalis bulk RNA, T1 Tempus bulk RNA and PacBio,
and T3 ONT tagged/deduplicated plus two Illumina single-cell products.
The six original T1/T2 queries are retained with their checksums and unchanged
default allele counts. This is not every public RNA library: native GRCh37,
unindexed products and further reprocessings remain outside scope.

## Sample, technology and assembly are separate axes

Every new protein comparison identifies its source product, then shows
assembly off/on and the recovered sequence/frame alternatives alongside
Varcode. Varcode is reference transcript plus the nominated edit, not measured
RNA or a unique transcript assignment. Missing reference predictions stay
visible. No callable reads, no alternate support, and no translated protein
despite alternate reads are different outcomes.

The protein result cap is explicitly removed **for these comparisons only**.
The default protein cap is unchanged. The old `mode.protein` top-result field
remains; `mode.proteins` retains every returned result, with all contributing
cDNA/frame contexts. These are recovered, reference-supported alternatives
under the recorded filters, not every biologically possible ORF or splice path.

Many uncapped results are shorter versions of another recovered protein.
Only exact, mutation-offset-aligned subsequences with the same strand, phase,
transcript assignments and frameshift status are grouped. A short ambiguous
sequence can belong to several alternatives. No longer protein is synthesized,
no template counts are summed, and every original rank remains in evidence
JSON. Pagination never discards an alternative; Varcode groups repeat alongside
RNA. Magenta marks differences where both tracks contain sequence, including
outside the focal mutation. Missing context is not a mismatch; this is an
edit-relative display, not a protein-homology alignment.

The earlier selected-read examples also retain all returned alternatives and
receive source labels. Their scope remains selected fixtures, not full-library
counts. The old Illumina protocol product used for DIAPH1 and one DYNC1H1
example has a conflicting path/date history; its timepoint is explicitly
unresolved rather than guessed from the filename.

## Additional small-indel support

Default alternate / reference RG-QNAME template counts in new products:

| Product | GTF3C5 | RNF213 | GLIS3 | KTN1 |
| --- | ---: | ---: | ---: | ---: |
| T0 BostonGene Illumina | 7 / 22 | 0 / 175 | 0 / 57 | 0 / 1048 |
| T0 Personalis Illumina | 0 / 227 | 0 / 274 | 0 / 54 | 0 / 1152 |
| T1 Tempus Illumina | 0 / 36 | 0 / 60 | 0 / 10 | 0 / 200 |
| T1 PacBio | 19 / 48 | 2 / 34 | 0 / 0 | 0 / 220 |
| T3 ONT dedup | 6 / 583 | 1 / 1047 | 0 / 1 | 0 / 1127 |
| T3 Illumina scRNA | 1 / 44 | 1 / 251 | 0 / 0 | 0 / 154 |
| T3 CD45neg Illumina | 0 / 17 | 1 / 38 | 0 / 1 | 0 / 82 |

Other/conflicting observations remain in the fixture; denominators above are
not total coverage or VAF estimates. Many PacBio input alignments lack QUAL;
Isovar 1.18.0 retains these reads with unknown base quality and unchanged
alignment/sequence filters. `--require-base-qualities` (API:
`ReadCollector(use_reads_without_base_qualities=False)`) restores the previous
strict behavior; `strict_quality_counts` retains that comparator in JSON.
No MAPQ, consensus accuracy or made-up constant is substituted for base Phred.
Known qualities remain unchanged. If overlapping mates disagree and either
base's quality is unknown, their evidence stays separate rather than guessing
which base wins. This default is centralized for CLI and API.

GTF3C5 reconstructs 49 aa in T0 BostonGene and T1 PacBio, and 32 aa in T3 ONT
in both modes. RNF213 reconstructs 41 aa in T1 PacBio in both modes. The new
PacBio contexts pass independent translation checks and agree with the
single-edit baseline. Previously all 19 GTF3C5 and two RNF213 alternate
template IDs were excluded solely because QUAL was absent.
The single alternate T3 Illumina observations do not pass the current
reconstruction support floor. Every returned protein in both modes passes
independent frame, cDNA translation and mutation-interval checks. All top
small-indel proteins match the single-edit expectation over recovered context.
Some **lower-ranked GTF3C5 ONT alternatives differ from that expectation**.
They remain visible, not relabelled as confirmed isoforms: ONT errors,
nonfocal variation and processing/molecular dependence are not adjudicated by
successful translation. RNF213, GLIS3 and KTN1's recovered alternatives agree
with the single-edit baseline.

## Rearrangements and larger deletion footprints

GABBR1--SLC29A1 has one additional complete tagged-ONT T3 path. Both actual
alignments pass MAPQ20 and the same-segment path check, but **none of its
20/30/45/60-nt-flank windows passes the extractor's Q10 requirements**.
This is path evidence, not a quality-validated new fusion sequence. No coding
frame is assigned. OTUD7A--FMN1 adds no complete paths in the new products;
the prior eight T2 paths/six CB-UMI labels remain unchanged.

A closer inspection separates two T3 GABBR1 limitations. The actual donor
alignment retains chr6:29612920, 14 bases from the catalogue's 29612906; the
catalogue endpoint itself is not mapped. Even using the observed endpoint,
the 20-nt-flank junction window contains Q3 bases (median Q27). The full read
has variable, valid Q3--Q50 scores (median Q40); its sequence and qualities
are identical in the tagged and deduplicated products. This is not evidence
of a garbled encoding. A diagnostic with synthetic high qualities still
fails at the catalogue endpoint and passes at the observed endpoint; it is
not accepted RNA evidence and no scores are changed in any fixture.

Across the bounded multi-locus RNA slices, the PacBio product has missing
QUAL on 720 of 727 primary MAPQ20 records. Its header identifies Iso-Seq 4.0
`groupdedup` followed by pbmm2. Missing QUAL after groupdedup is a
[previously reported workflow limitation](https://github.com/PacificBiosciences/pbbioconda/issues/694),
not proof of an aligner bug or poor underlying HiFi reads. Recovering earlier
quality-bearing reads can add calibrated base confidence; retaining the
current consensus sequence without fabricating Phred avoids losing it. Other
queried products have QUAL on every primary MAPQ20 record; no sequence/QUAL
length mismatch was found. Illumina's small discrete score sets are consistent
with [quality-score binning](https://emea.support.illumina.com/content/dam/illumina-support/documents/documentation/system_documentation/novaseq/1000000019358_18_novaseq-6000-system-guide.pdf),
not by themselves corruption. These counts describe queried slices, not the
entire source BAMs.

The original missing-QUAL PacBio records retain `ic`/`is` consensus/read-count
tags, not calibrated per-base confidence; those counts are not independent
molecular support. The existing MT-ND5 corpus case also regains a translated,
independently validated protein. KTN1's PacBio reference count becomes 220;
the NTF3 compound case regains one reference and one alternate observation
but remains below the reconstruction support floor. Original BAMs are unchanged.

The [soft-clip on/off audit](SOFT_CLIPS.md) keeps missing-QUAL policy fixed and
does not change the clipping default. Extra unaligned ends did not rescue a
protein in the four-indel panel and often reduced supported context.

AFF3 and KEAP1 retain ordinary annotated splice skips across their intronic
DNA intervals in additional products. These normal mature-RNA paths cannot
distinguish the intronic deletion from undeleted DNA. All three large events
have no matching RNA D/N operation or observed supplementary adjacency at
their originally listed endpoints. DLG5's expanded endpoint and sequence
checks are below. RNA within an interval is retained-region expression, not
coverage-normalized expression or evidence against a subclonal DNA event.

## DLG5: resolve DNA first, leave mutant RNA unresolved

The [DRAGEN source calls](https://osteosarc.com/dragen/sv/) contain more than
the 79,531-bp deletion label: `SVINSSEQ=GAAATGATGC`, a 10-nt breakpoint insert.
The [DRAGEN specification](https://help.dragen.illumina.com/product-guides/dragen-v4.3/dragen-dna-pipeline/sv-calling)
describes this field separately from the deletion and defines CONTIG as the
assembled sequence. Three source records (T1, T2, organoid) and their original
headers are pinned; T0 has no corresponding call in the queried VCF.

Purple/ESVEE instead records an explicit BND at chr10:77850914 with
`GCTTCTCTGAAATGATGCTTCTCCA[chr10:77930461[`. Removing its reference anchor
gives 24 inserted bases between retained boundaries 77850914 and 77930460.
That full sequence with GRCh38 flanks occurs in all three DRAGEN CONTIGs and
in original tumor DNA reads. Concatenating the nominal DRAGEN DEL/insert
fields onto reference flanks does **not** reproduce the complete haplotype:
nearby sequence changes matter. This is not a claim that the focal DRAGEN
call is invalid, nor that its endpoints can silently be rewritten.

Using the explicit BND plus exact 20-nt reference flanks, Q20 primary/MAPQ20
sequence support is 2 T1 DNA templates, 4 T2 and 4 organoid; none in the queried
T0 tumor or two normal products. Actual compatible supplementary paths number
3, 4 and 4 respectively, with CIGAR boundaries 77850921/77930460. These are
overlapping evidence categories and not interchangeable with caller PR/SR
counts. The original DRAGEN calls report tumor alternate PR/SR of 3/7 (T1),
4/5 (T2), 2/5 (organoid), and zero alternate PR/SR in their matched normal.
The BND source, assembled contig, alignment and strict exact-window counts
answer different questions; all representations are preserved.

No exact junction signature is recovered in any of the 14 queried RNA
products, even in the explicitly unqualified sequence-only check. No matching
D/N or observed supplementary path is found at the original, explicit-BND or
observed-DNA-CIGAR endpoints. This cannot exclude a spliced RNA product that
removes the intronic junction, low expression or subclonality.

The event covers 864 coding bases and the annotated start of complete
DLG5-001 (Ensembl 87). Its incomplete DLG5-002 model has no annotated start.
Q20 RNA retaining the normal DLG5-001 start is observed in T1 bulk (6 templates),
T2 bulk (4), T0 Personalis (90), T1 Tempus (6), T3 scRNA (24) and T3 CD45neg (22).
Six PacBio templates span this start without adequate quality evidence.
Retained-start RNA and mutant DNA can coexist: neither an expressed mutant
CDS nor its absence follows. No guessed protein or junction peptide is emitted.

Reads can therefore be fetched and their splice footprints inspected even
when Varcode has no concrete protein prediction. Ordinary/novel splice paths
are **candidate consequences**, not mutation-assigned reads, until direct or
safely phased sequence evidence links them to the mutant allele. Annotation
and presence alone cannot supply that link.

## DNA-nominated intergenic partners and retained clips

Three literal T1 Purple BND records nominate two WIPF2 breakends
(chr17:40221092 to chr9:42953930; chr17:40221238 to chr6:9473195)
and TMEM63B (chr6:44154877 to chr12:133264867). Ensembl 87 has no gene at
the partner positions. The WIPF2 records share an assembly identifier and
assembly links: they are components of a complex event, not two independently
established fusions. Varcode 9.2.3 classifies coding overlaps as
`TranslocationToIntergenic` without a concrete protein. That effect name is
broader than a literal intergenic partner; here the partner annotation was
checked independently. Both BND sides must be queried, since the mate of an
intergenic-anchored record can be genic.

DNA alone supplies useful search constraints: retained sides/orientation,
inserted junction sequence when present, affected transcripts/exons and
whether annotated coding starts/stops are retained. These nominate RNA
queries and candidate splice paths, not an expressed CDS. An ordinary splice
can remove the DNA junction, so footprint linkage remains necessary.

Original T1 tagged ONT, bulk Illumina and PacBio regional BAMs were acquired
for the complete two genes and +/-5-kb partner windows, with matching GRCh38
reference windows. An explicit screen retained all >=25-nt primary/MAPQ20
soft clips, including missing-QUAL PacBio clips. It searched both orientations
for diverse exact 25mers unique **within that partner window**, not genome-wide.
No observed compatible supplementary path or SA nomination to these partner
windows was found. TMEM63B has two ONT and three Illumina clips within 100 bp
of the DNA breakend, but none matches a screened partner seed.

WIPF2 has no >=25-nt clip within 100 bp of either nominated breakend.
Across its whole gene, three ONT clips have 27-base matches to the chr9 window;
six ONT and one missing-QUAL PacBio clip have 25--27-base matches to chr6.
Several identical clipped reads match both partner windows, and clipping
boundaries lie 1.2--49.8 kb from the DNA breakends. These are weak candidate
matches, not unique partner alignments or mutation-assigned RNA. No frame or
protein is inferred. The bounded exact-seed screen can miss error-bearing or
spliced partners; absence here is not biological absence.

The dated output retains raw BAMs, indices, exact VCF records, reference
snapshots, every clip summary/hit SAM and the screen script. A justified next
step is splice-aware, competing-placement realignment of nominated clips,
not globally treating every soft clip as coding context.

## Provenance and reproduction

ONT tagged/deduplicated files share processing families and are not added.
T3 scRNA/CD45neg and long-read products are not established independent
molecular replicates. Labels, UMIs and different filenames do not prove
independence or malignant-cell origin. Source metadata is linked through the
hash-pinned NR2F2 library inventory; acquired BAM/index hashes, URLs, query
regions and dates are retained. Query windows are bounded; no whole BAMs were
downloaded. Failed sandbox/network attempts are not biological negatives.

```sh
python -m tests.data.osteosarc.figure_comparisons.extended_footprints \
  --acquire /new/rna --baseline-inputs /previous/run/inputs
python -m tests.data.osteosarc.figure_comparisons.extended_footprints \
  --inputs /new/rna --output /new/rna-audit
python -m tests.data.osteosarc.figure_comparisons.dlg5 --acquire /new/dna
python -m tests.data.osteosarc.figure_comparisons.dlg5 --context /new/context
python -m tests.data.osteosarc.figure_comparisons.dlg5 --audit /new/dlg5-audit \
  --dna /new/dna --rna /new/rna --context-inputs /new/context --footprints /new/rna-audit
python -m examples.osteosarc_context_figures --output-dir figures/osteosarc
isovar plot --bam sample.bam --variant 9 133057893 GGAGGAGGAGGAA G \
  --genome GRCh38 --compare-assembly --all-proteins --sample-label 'T1 / ONT'
```

The generated audit can be pinned with `extended_footprints --pin <audit>
--pin-name extended-footprints` (or `dlg5`); existing pins are never silently
overwritten. Tests recount DLG5 start/signature support and supplementary paths
from original pinned SAMs, and check every protein/frame's validation and
display coverage. Complete indel/large-interval raw recounts still require the
retained regional BAMs or reacquisition, not just their summary JSON.
