# Nine additional osteosarc RNA footprints

Audited September 17, 2026 against original public GRCh38 RNA alignments and
Ensembl 87 annotation. This is a bounded T1/T2 analysis, not a claim about all
available timepoints or libraries. Sources: [indel catalogue](https://osteosarc.com/oncoanalyser/snv/),
[SV catalogue](https://osteosarc.com/dragen/sv/), and
[fusion catalogue](https://osteosarc.com/fusions/).

## Explicit small indels

Default Isovar alternate-template counts (RG/QNAME identity):

| Gene | GRCh38 VCF allele | T1 ONT | T2 ONT | T1 short | T2 short |
| --- | --- | ---: | ---: | ---: | ---: |
| GTF3C5 | chr9:133057893 GGAGGAGGAGGAA>G | 117 | 128 | 7 | 12 |
| RNF213 | chr17:80327830 ATAC>A | 5 | 14 | 42 | 11 |
| GLIS3 | chr9:3856149 CTGATGTGG>C | 0 | 0 | 9 | 0 |
| KTN1 | chr14:55627965 G>GTT | 0 | 10 | 0 | 2 |

Protein panels deliberately exclude secondary alignments; defaults are not
changed. This yields GTF3C5 short-RNA alternate counts 6/10, and KTN1 T2 short
count 3. For KTN1, one primary-supporting template becomes conflicting/other
when its alternate placement is considered. Default and primary-only counts
are both retained, not silently substituted. All other listed alt counts agree.

Both assembly modes pass independent frame, RNA translation, edited-reference
and mutation-interval checks for every reconstructed top protein. All four
agree with the single-edit reference prediction over recovered context:
GTF3C5 removes four glutamates, RNF213 removes leucine, GLIS3 and KTN1 shift
frame. Residue numbers vary by transcript/repeat normalization; literal
alleles and full sequences are retained.

| Gene/product | Assembly on / off, aa |
| --- | ---: |
| GTF3C5 T1 ONT; T1 short | 34/34; 49/47 |
| GTF3C5 T2 ONT; T2 short | 32/32; 49/49 |
| RNF213 T1 ONT; all other assessed products | 29/29; 49/49 |
| GLIS3 T1 short | 46/46 |
| KTN1 T2 ONT; T2 short | 29/29; 18/18 |

GLIS3 has only 2/1 reference-supporting ONT templates and no callable T2 short
template: no protein is not absence of the allele. KTN1's ONT reconstruction
extends farther upstream than short RNA, but differing support/library
composition prevents attributing that difference solely to read length.

## Rearrangement RNA

- **GABBR1--SLC29A1**, chr6:29612906 / chr6:44218908: 2 T1 and 3 T2
  complete tagged-ONT paths. Each sample supports one exact 120-nt Q10 window
  with a direct junction and no inserted bases. Neither side matches an
  exact collinear annotated transcript at this junction: frame unresolved.
- **OTUD7A--FMN1**, chr15:31743670 / chr15:33043479: 8 T2 paths, six distinct
  CB/UMI labels. Four read/label observations support the selected exact
  122-nt Q10 window, including a two-base `AG` insert in the displayed -/+
  orientation. FMN1 coding annotation is on the opposite strand; the OTUD7A
  breakpoint is outside its complete coding isoform. No coding frame assigned.

Both actual pieces must pass MAPQ20 and reciprocal SA/path compatibility.
No high-MAPQ complete path is observed in the two queried short-RNA products.
Deduplicated ONT predecessors retain SA declarations but lack complete
matching partner records: their zero complete-path count is **not** absence
of fusion RNA. Tagged and deduplicated records must never be added. Neither
QNAME counts nor CB/UMI labels establish independent molecules or malignant
cell identity. No long-read-only biological claim or fusion protein is made.

## Larger DNA deletions

- **DLG5**, chr10:77850921-77930452 (79,531 bp): overlaps 864 coding bases
  of DLG5-001 and 639 of DLG5-002. It includes the annotated start of complete
  DLG5-001; DLG5-002 is incomplete and has no annotated start codon.
  This is not a simple internal in-frame deletion merely because those lengths
  are divisible by three. No matching RNA D/N or observed supplementary
  adjacency was recovered. Aligned RNA inside the interval: 100/81 ONT and
  314/454 short-RNA template IDs, T1/T2. This supports retained-region
  expression, not exclusion of a subclonal event or expression quantitation.
- **AFF3**, chr2:99668095-99668223 (128 bp): no annotated exon overlap.
  The ordinary intron 99649666-99672537 spans the DNA interval: 6/16 ONT and
  17/2 short-RNA templates, T1/T2. Mature RNA taking this splice path cannot
  distinguish the intronic deletion from the undeleted allele.
- **KEAP1**, chr19:10492683-10492894 (211 bp): no coding-exon overlap; 34
  bases overlap a processed-transcript exon (KEAP1-007). The ordinary intron
  10492262-10499394 spans the interval: 473/535 ONT and 71/23 short-RNA
  templates. Again, this ordinary splice is not evidence for the DNA deletion.

No matching D/N event or observed supplementary adjacency at the listed ends
(within 3 bp) occurs for these three deletions in the acquired products.
Coverage is aligned M/=/X bases, never CIGAR N/D reference span. Novel splice
paths remain in the JSON but are not attributed causally to the DNA event.
Counts are primary/QC-passing, nonduplicate, MAPQ >=20 (255 retained);
structural path checks also assess tagged supplementary records. No RNA allele
is fabricated from symbolic SV placeholders or partner labels.

## Reproduce and inspect

```sh
python -m tests.data.osteosarc.figure_comparisons.footprints --acquire /new/inputs
python -m tests.data.osteosarc.figure_comparisons.footprints --inputs /new/inputs --output /new/audit
python -m tests.data.fusions.build_osteosarc --alignments /new/inputs --output /new/fusions
python -m examples.osteosarc_context_figures --output-dir figures/osteosarc
```

Acquisition downloads only indexed regions, never whole BAMs. Complete source
URLs, regions, commands, acquisition dates and regional BAM/index hashes are
in receipts. Six input BAMs are retained with the delivered figure run, not
duplicated into git. The compact hash-pinned `rna-footprints.json.gz` contains
all summaries, model intervals, independent checks and original SAMs for
new chimeric paths. Indel/large-deletion counts require the retained input
BAMs or reacquisition; they are not falsely described as offline raw-read
regressions. Tests independently recount chimeric paths/windows and verify
translation from pinned cDNA, with synthetic D/N and missing-partner controls.

Individual panels are white-background 600-dpi PNG and vector SVG; each
example has its own PDF, with the full bookmarked `isovar-all-figures.pdf`
preserving earlier examples. Protein order remains no assembly, assembly,
Varcode. Transcript identifiers include names. Counts/limitations sit outside
the core tracks. No new public Isovar API or default is introduced.
