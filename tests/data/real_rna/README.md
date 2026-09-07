# Public RNA regression fixtures

## Scope and verification plan

Preserve small, unmodified public RNA alignment sets across platforms and
aligners; compare allele extraction with an independent CIGAR ledger and DNA
benchmark; expose quality-policy assumptions before implementing weighting.
No production quality cutoff, weighting, or ranking change is shipped here.
The earlier Q20 implementation draft is deferred, not validated by this PR.

These are regression examples, **not** a representative accuracy study or a
clinical validation set. NA12878 variants are germline DNA variants, not
confirmed somatic neoantigens. PacBio K562 reads are alignment/quality edge
cases; no small-variant truth is asserted for that sample. The fixtures do not
validate translations, peptide ranking, or current sequencing chemistries.

## Sources and selection

Accessed 2026-09-07. All coordinates are GRCh38. SAM headers, sequences,
qualities, names, flags, CIGARs and auxiliary tags are retained. We only select
records and serialize/compress them; there is no realignment, base editing,
quality replacement, downcoding, or added MD/NH tag. BAM-to-SAM conversion can
change binary encoding, so hashes pin the resulting SAM text, not source BAM
bytes. `manifest.json` also identifies every individual record by SHA-256.

| Fixture | Verified platform / source | Aligner in original header | Selection |
| --- | --- | --- | --- |
| `illumina_star.sam.gz` | Illumina HiSeq 2000, paired-end NA12878 RNA, SRR1258218 | STAR 2.7.10b | 122 records: four chr4 loci, available mates/alternative alignments, three multimapped templates, one overlapping pair |
| `ont_minimap2.sam.gz` | Oxford Nanopore MinION direct RNA, NA12878; Albacore 2.1-era data | minimap2 2.5-r572 | All 392 alignments overlapping chr6:43767094–43788458 |
| `pacbio_minimap2.sam.gz` | PacBio Sequel / Sequel II K562 Iso-Seq, ENCODE via JBrowse | minimap2 2.15-r905 | 28 of 579 regional records: first four per (FLAG, set of QUAL scores) group |

The PacBio selection is deliberately stratified to retain both orientations,
supplementary alignments and quality conventions. Do not infer frequencies
from it. Regional subsets do not necessarily include every mate or distant
secondary/supplementary alignment; NH and SA remain as originally reported.

### Illumina

- [ENA SRR1258218](https://www.ebi.ac.uk/ena/browser/view/SRR1258218), sample
  SAMN02731489 / GSM1372331, study PRJNA245078. ENA identifies the sample as
  `NA12878_RNASeq` and instrument as `Illumina HiSeq 2000`.
- [Kent Riemondy's raerdata preparation recipe](https://github.com/rnabioco/raerdata/blob/devel/inst/scripts/NA12878-data-processing.Rmd)
  and [Bioconductor dataset documentation](https://www.bioconductor.org/packages/devel/data/experiment/vignettes/raerdata/inst/doc/raerdata.html).
- RNA BAM: [ExperimentHub EH8466](https://experimenthub.bioconductor.org/fetch/8466);
  index: [EH8467](https://experimenthub.bioconductor.org/fetch/8467).
  BAM MD5: `dc15bcd56e25e8ad09937fe64594438a`, 2,262,239 bytes.
  Resolved mirror:
  `https://mghp.osn.xsede.org/bir190004-bucket01/ExperimentHub/raerdata/NA12878/1.0.0/NA12878.rnaseq.sub.bam`.
- The published derivative is already restricted to chr4:1–1000000 and has
  Picard-marked duplicate records removed with `samtools view -F 1024`. It is
  not suitable for testing duplicate prevalence. Original @PG paths say
  `NA12877`; those paths are preserved, but sample attribution follows the
  actual SRR1258218 archive metadata and published dataset description.

### Nanopore

- Workman et al., *Nanopore native RNA sequencing of a human poly(A)
  transcriptome*, Nature Methods (2019), [doi:10.1038/s41592-019-0617-2](https://doi.org/10.1038/s41592-019-0617-2).
- [Consortium RNA data README](https://github.com/nanopore-wgs-consortium/NA12878/blob/master/RNA.md).
- Original indexed BAM:
  `https://s3.amazonaws.com/nanopore-human-wgs/rna/bamFiles/NA12878-DirectRNA.pass.dedup.NoU.fastq.hg38.minimap2.sorted.bam`.
  Index is the same URL plus `.bai`. The upstream filename indicates pass,
  dedup and NoU preprocessing. Original command is `minimap2 -ax splice -uf
  -k14`; no MD or NH tags are supplied in this region. This old basecaller is
  not a model of present-day nanopore performance.

### PacBio

- [ENCODE ENCSR589FUJ](https://www.encodeproject.org/experiments/ENCSR589FUJ/)
  and [ENCSR983KDL](https://www.encodeproject.org/experiments/ENCSR983KDL/).
  Instrument metadata is on input files ENCFF694INI / ENCFF763VZC (Sequel)
  and ENCFF429VVB / ENCFF634YSN (Sequel II).
- [JBrowse K562 example](https://jbrowse.org/jb2/docs/tutorials/k562_fusions/)
  exposes the indexed BAM
  `https://jbrowse.org/demos/cancer_sv/K562_isoseq.bam` (`.bai` index).
  Its @PG chain records sorting/merging ENCODE unfiltered genomic BAMs
  ENCFF433YKW, ENCFF092NLB, ENCFF515YRZ and ENCFF475XQX. We extract
  chr22:23286000–23293000. These are not TranscriptClean-corrected alignments.
- No per-record PG/RG tag associates a read with one of the four input files;
  we therefore label the fixture as the pool, not an individual instrument run.

### Independent DNA benchmark

[GIAB HG001 / NA12878 v4.2.1 GRCh38](https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/NA12878_HG001/NISTv4.2.1/GRCh38/),
`HG001_GRCh38_1_22_v4.2.1_benchmark.vcf.gz` and corresponding `.bed`.
`HG001.selected.vcf` contains eight original records and the original VCF
header. `HG001.selected.bed` contains original benchmark intervals covering
every selected allele, checked by a test.

| Chromosome | Original VCF POS / REF / ALT |
| --- | --- |
| chr4 | 337869 T C; 337877 G C; 337904 G GC; 768474 GTT G |
| chr6 | 43770613 C G; 43780233 G GGT; 43785475 A G; 43785588 G A |

Matched cell-line DNA supports the locus labels; it does not make every RNA
read correct or establish a tumor mutation. RNA editing, culture differences,
mapping errors, and basecalling errors are still possible.

## What these data establish

- The selected STAR records use MAPQ 255 exactly when NH=1, while multimapped
  records have MAPQ 0 or 3. [SAM](https://samtools.github.io/hts-specs/SAMv1.pdf)
  defines 255 as unavailable, not Phred 255; [STAR's own parameters](https://github.com/alexdobin/STAR/blob/master/source/parametersDefault)
  specify its unique-mapping convention. NH is an alignment count, not a
  calibrated probability or an independent count of RNA molecules.
- At chr6:43785588, 106 primary nanopore alignments have the alternate **base**
  A: 63 have base quality >=10, 4 >=20, and none >=30. These are base-only
  observations, not all pure SNP allele support: adjacent insertions must be
  retained as compound alleles. A universal Q20 cutoff would discard most
  alternate-base observations in this example. This does not establish which
  rejected reads are accurate.
- At chr4:337869, two of eight primary alternate observations have quality Q2.
  High mapping confidence alone therefore does not guarantee high base quality.
- All 579 PacBio regional records have a constant quality within each read,
  either Q20 or Q40. These supplied values cannot discriminate individual
  bases within a read. The fixtures do not establish whether the scores are
  calibrated consensus confidence or pipeline placeholders.
- STAR's real repeat-shifted insertion has different qualities at the original
  and normalized query positions. Tests retain both values; future quality
  filtering must not silently substitute one for the other.
- **Known bug [#215](https://github.com/openvax/isovar/issues/215):** CIGAR N can
  be returned as a deletion. Three strict expected-failure tests preserve the
  reproducer across all platforms. This PR does not fix or conceal that bug.

These observations inform [#8](https://github.com/openvax/isovar/issues/8),
but do not implement quality weighting or choose new default thresholds.
Future work needs explicit unknown/uninformative-quality behavior, an
aligner-aware mapping policy, fragment-level counting, and sensitivity checks
that do not turn a quality score into a claimed variant posterior.

## Offline tests and reproducibility

Run `pytest -q tests/test_real_rna.py`. Tests use the checked-in reads only and
create disposable BAMs/indexes. No network, reference genome, gene annotation,
STAR, minimap2 or external samtools executable is required for these tests.
The CIGAR ledger generator does not import isovar. It distinguishes D from N,
retains adjacent insertions, and records per-read query positions and qualities;
normalization is checked separately against an explicitly inspected real read.

To independently download and regenerate (requires network, pysam and samtools;
use a new source directory):

```sh
python tests/data/real_rna/fetch_sources.py /tmp/isovar-public-rna-sources
python tests/data/real_rna/rebuild.py /tmp/isovar-public-rna-sources /tmp/isovar-public-rna-rebuilt
diff -r tests/data/real_rna /tmp/isovar-public-rna-rebuilt
```

README/scripts/license are not generated, so the last command reports those
as present only in this directory. Generated data/manifest/VCF/BED should match.
The manual fetch downloads a 2.3-MB STAR BAM, regional ranges of indexed long-read
BAMs, regional VCF records/index and a 15.5-MB benchmark BED. It does not download
the full long-read BAMs. Upstream URLs may change; do not bless a changed hash
without inspecting the new source and ledger. Deterministic gzip uses mtime 0;
the integrity tests use uncompressed hashes to avoid compressor-version issues.

## Reuse and attribution

These are third-party public research data, not newly generated isovar data.
The nanopore consortium explicitly releases these RNA data under **CC-BY**;
retain attribution to Workman et al. and the consortium and identify our
regional subsetting. See its linked RNA README. ENCODE permits reuse with
attribution under its [data policy](https://www.encodeproject.org/help/faq/);
cite the experiments/accessions above and JBrowse for the indexed derivative.
ENA's [public data policy](https://www.ebi.ac.uk/ena/browser/about/policies) and
the linked GIAB release describe the other public sources. raerdata is MIT
licensed, copyright 2023 RNA Bioscience Initiative; its notice is retained in
`RAERDATA-LICENSE.txt`. No controlled-access patient data were used.
