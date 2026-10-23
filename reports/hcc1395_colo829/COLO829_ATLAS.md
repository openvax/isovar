# COLO829 whole-run original-RNA atlas

This report joins the matched-normal SMaHT COLO829/COLO829BL call set to the
complete public Illumina RNA run SRR8615617. The RNA comes from a separate
COLO829 culture, so an RNA reference allele does not negate a SMaHT DNA
somatic call. TAGET polished representatives and full-length cluster counts
remain in the separate tables documented in [README.md](README.md); those
counts are not original allele-spanning read support.

## Inputs and processing

- SMaHT source inventory: 44,005 SNVs and 2,059 indels. The source VCF files
  and SHA-256 values are listed in [README.md](README.md#input-sha-256).
- Separate SV benchmarks: 54 Paulin Paz 2024 v2 VCF records and 134
  Valle-Inclán 2022 lifted VCF records. A VCF breakend record is a record,
  not necessarily an independent event. The two benchmarks are not combined.
- [SRR8615617](https://www.ncbi.nlm.nih.gov/sra/SRR8615617) is paired
  Illumina cDNA RNA-Seq, BioSample SAMN10988208, study PRJNA523380. The
  [ENA file report](https://www.ebi.ac.uk/ena/portal/api/filereport?accession=SRR8615617&result=read_run&fields=run_accession,sample_accession,study_accession,library_strategy,library_source,library_layout,read_count,base_count,fastq_bytes,fastq_md5&format=tsv)
  lists 115,269,780 paired spots (230,539,560 individual 101-base reads) and
  23,284,495,560 bases. SRA Toolkit `sra-stat --quick` independently confirms
  the spot and base totals.
- The full 12,977,933,218-byte [NCBI normalized SRA archive](https://sra-pub-run-odp.s3.amazonaws.com/sra/SRR8615617/SRR8615617)
  has SHA-256 `9d89bb4f5a74fba7a8ba9495fd44046c35f4abb60d10c7bd9261f07569520b4f`
  and was validated with SRA Toolkit 3.4.1. `fastq-dump -Z --split-spot` streams
  every biological read pair to minimap2 2.26 (`-ax sr --frag=yes`) against
  the NCBI GRCh38 no-alt analysis set. The FASTA gzip SHA-256 is
  `fb4243ebb014caf27111f24dd62b7ce42160f28581da6f8fcd6cba5977778d02`.
  The complete run is aligned; an indexed BAM retains alignments overlapping
  any source small-variant allele or a 500-base window at either SV end.
  The merged target BED has 45,364 intervals covering 169,087 reference
  bases. This is output selection after whole-run alignment, not a read-ID
  preselection from the archive's older GRCh37 alignment. The `sr` preset is
  genomic rather than splice-aware; junction-crossing reads can be lost or
  placed incompletely. Zero support therefore applies to this alignment
  method and cannot establish absence of an expressed allele.
- The indexed 10 MB selected-site BAM and exact SV source VCFs are committed in
  `tests/data/colo829_whole_run/`. Their SHA-256 values and scope are listed
  in the [fixture README](../../tests/data/colo829_whole_run/README.md).
  `fastq-dump` reported 115,269,780 spots read and written; summing minimap2's
  466 batch logs yields exactly 230,539,560 individual reads. The selected-site
  BAM has 129,534 alignment records (125,980 primary) and passes
  `samtools quickcheck`.
- Isovar 1.59.6 with Ensembl 87 counts allele-read objects after mate
  handling; `*_fragments` counts distinct template identities. These are not
  UMI-confirmed independent molecules. Translation uses Isovar defaults and
  retains all ranked local contexts through 750 ALT allele-read observations.
  Above 750, it evaluates all ALT reads as read-spanning local contexts at
  one context length, reconciles compatible-read support, trims to Isovar's
  default coverage floor, and retains all ranked protein contexts. It does
  not overlap-assemble different read contexts at these loci. The threshold
  chooses a computational method, not a biological status; no ALT-supported
  locus is excluded from translation solely due to read depth. Full allele
  and fragment counts remain in the call table.

## Whole-run findings

All 46,064 SMaHT small-variant source rows have an RNA status. Among 44,005
SNVs, 1,865 have at least one observed ALT allele; among 2,059 indels, 68 do.
These are allele observations in this separately cultured RNA sample, not
validated full-length mutant isoforms. The detailed status breakdown is:

| RNA status | SNVs | Indels | Total |
|---|---:|---:|---:|
| `alternate_observed` | 1,865 | 68 | 1,933 |
| `reference_only` | 1,515 | 122 | 1,637 |
| `reference_and_other_only` | 33 | 3 | 36 |
| `other_allele_only` | 23 | 4 | 27 |
| `alignment_without_allele_call` | 149 | 17 | 166 |
| `no_alignment_overlap` | 40,420 | 1,845 | 42,265 |

The 1,933 ALT-observed calls have 35,949 ALT allele-read observations in
total, represented by read IDs in the witness table. Local translation
produced 12,218 ranked protein contexts for 124 calls, including four indels.
The other 1,809 ALT-observed calls had no local translation under the stated
method. Multiple contexts from one variant are not resolved full-length
isoforms. Three calls used the high-depth read-spanning method: at
`chr9:21971154:CCG>C`, all 6,984 ALT allele reads were evaluated, yielding
1,416 local contexts, including a CDKN2A frameshift context with 5,338
supporting reads; `chr2:105361307:G>GCCCC` (847 ALT reads) yielded 146
contexts; `chr7:101210439:G>A` (1,300 ALT reads) yielded no local coding
context. No local translation does not imply absence of a mutant protein.

The 188 source SV VCF records have the following screen outcomes:

| RNA status | VCF records |
|---|---:|
| `no_breakend_window_alignment` | 78 |
| `one_end_covered_no_split_candidate` | 36 |
| `both_ends_covered_no_split_candidate` | 57 |
| `cross_end_pair_context_unresolved` | 11 |
| `unresolved_partner_coordinate` | 6 |

No record had an observed same-read split-alignment candidate in the selected
BAM. The 11 mate-context records contain 26 record-level candidate counts;
the two VCF benchmarks and paired breakend records can count the same SRA
read more than once. For example, Paulin Paz v2
`COLO829_SOMATIC_SV3_DEL` has 427 and 956 read names in its two windows and
six opposite-mate pair names, but zero observed split-alignment candidates.
Those pairs do not establish its deletion junction or mutant protein.

The [selected junction audit](sv_junction_audit.tsv) finds three additional
SV3 split-junction fragments only when Isovar accepts alignments down to
MAPQ 1. Their linking pieces have MAPQ 1-6, so none survives MAPQ 20.
SV42 and SV49, the two exon-involved v2 deletions with mate candidates,
have no junction-spanning read at either threshold in this selected-site BAM.
The corresponding earlier truth-set records reuse the same mate names.
The original paired RNA was aligned in genomic short-read mode; its inability
to resolve a junction here is not proof that the RNA isoform is absent.

## Outputs

- `colo829_atlas_calls.tsv.gz`: one row for each of the 46,064 SMaHT small
  variants. It preserves the DNA source key and filter, raw alignment overlap
  for otherwise uncallable sites, reference/alternate/other read and fragment
  counts, ALT source-read count, genomic transcript overlaps, and explicit RNA
  and protein statuses. `chimeric_path_warning_reads` counts distinct read IDs
  for which Isovar rejected an inconsistent supplementary path; it does not
  imply those reads support an SV. `assembly_status=greedy_merge_skipped_complexity`
  identifies loci where Isovar's default 1,000-candidate guard bypassed its
  greedy sequence merge; the returned local contexts remain available but
  cannot be read as completed assemblies.
- `colo829_atlas_proteins.tsv.gz`: every ranked, locally translated mutant
  context from alternate-supporting reads, with zero-based mutation interval,
  supporting fragments/reads, and compatible reference transcript IDs. A
  local amino-acid context is not a full-length mutant isoform or protein
  expression observation. Multiple rows can share a residue but differ in
  available context; rank is Isovar's preference, not a calibrated posterior.
- `colo829_atlas_alt_witnesses.tsv.gz`: original SRA read IDs and source-read
  counts for each ALT observation. A read can appear for several nearby calls;
  cross-variant linkage requires its original sequence and path to be checked.
- `colo829_sv_rna_screen.tsv`: every VCF record in both separate SV benchmarks,
  preserving the VCF ALT and DNA support, with alignment and fragment counts
  in each 500-base breakend window. The
  `observed_split_alignment_candidates` count requires two mapped segments
  from the same read end, one supplementary, observed at the two windows.
  `cross_end_mate_pair_candidates` instead counts opposite primary mates at
  distinct windows; that is supporting context, not a junction-spanning read.
  This is a **candidate screen**. Its path and protein statuses remain
  unresolved because alignment proximity alone does not establish junction
  sequence, orientation, splice interpretation, reading frame, or a full ORF.

The small-variant `rna_status` values mean:

| Status | Meaning |
|---|---|
| `alternate_observed` | At least one Isovar ALT allele-read observation. |
| `reference_only` | Callable reference reads and no ALT or other allele. |
| `reference_and_other_only` | Reference and other bases, without the nominated ALT. |
| `other_allele_only` | Callable bases, but none match REF or nominated ALT. |
| `alignment_without_allele_call` | A genomic alignment overlaps, but Isovar cannot call an allele. |
| `no_alignment_overlap` | No retained alignment overlaps the reference span. |

`no_translation` means ALT reads were present but Isovar did not resolve a
default local coding context. It is not evidence against translation in the
cell. `not_attempted_no_alt` is distinct from a failed translation. The
`multiple_local_contexts` status retains ranked alternatives without claiming
that any one is a resolved full-length isoform. The
`genomic_overlapping_transcript_ids` column lists Ensembl models at the DNA
coordinate, not RNA-compatible isoforms; compatibility for translated
candidates appears in the protein table.

## Rebuild

The committed source inventory, two SV VCFs, and indexed selected-site BAM
allow offline regeneration of all four report tables with Ensembl 87 cached:

```sh
python reports/hcc1395_colo829/build_colo829_atlas.py
python -m reports.hcc1395_colo829.build_colo829_sv_screen
```

To recreate the selected-site BAM from the full public run, install SRA
Toolkit 3.4.1, minimap2 2.26, and samtools 1.21. Obtain the NCBI GRCh38
no-alt FASTA named in [README.md](README.md). Use a scratch directory with
space for the SRA archive, reference FASTA, and temporary alignment files.

```sh
curl -L -C - -o SRR8615617.sra \
  https://sra-pub-run-odp.s3.amazonaws.com/sra/SRR8615617/SRR8615617
python reports/hcc1395_colo829/build_atlas_regions.py \
  --output atlas-regions.bed
reports/hcc1395_colo829/align_full_colo829_rna.sh \
  GRCh38-no-alt.fa atlas-regions.bed SRR8615617.sra \
  /path/to/sratoolkit.3.4.1/bin/fastq-dump atlas-work
shasum -a 256 atlas-work/SRR8615617.all_sites.grch38.bam
```

The alignment script validates the archive before decoding it. Compare the
BAM hash with the fixture README before replacing any report table. Because NCBI
stores this SRA run as an aligned archive, SRA Toolkit may retrieve public
reference sequence while decoding it; the offline tests never do so.

## Unresolved interpretation

The SV screen does not yet reconstruct alternative breakpoint-spanning RNA
paths, compare normal splicing and other genomic matches, or infer an ORF.
That work belongs to [Isovar #305](https://github.com/openvax/isovar/issues/305)
and [#306](https://github.com/openvax/isovar/issues/306); #541 remains open for
the COLO829 application. HCC1395 whole-run RNA remains in
[#533](https://github.com/openvax/isovar/issues/533). Vaxrank follow-up is
tracked in [Vaxrank #658](https://github.com/openvax/vaxrank/issues/658).
