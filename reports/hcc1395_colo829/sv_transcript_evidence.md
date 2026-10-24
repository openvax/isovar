# Original RNA evidence for selected SV transcript paths

These are **bounded, breakpoint-compatible RNA path hypotheses**, not a
complete atlas of SV transcripts. The selected original-read fixtures cover
three HCC1395 and one COLO829 DNA events with direct junction observations.
They yield **7 retained paths** (6 HCC1395;
1 COLO829). Only **2**
have at least two exact local path witnesses. HCC1395 chr9 reaches
`path_limit`, so this is a retained-hypothesis count rather than a
count of biological isoforms.

**Evidence units:** direct fragments are distinct original RNA fragments
making the join, shared by paths. Local exact fragments each match one
reported window across the join with compatible genomic placements.
The 12-base-flank window is reported separately when it cannot
distinguish competing paths. A `unique_among_retained_paths` window
differs in sequence from other retained event paths; uniqueness is
limited by the path limit and does
not establish a full isoform. Exact matching can undercount noisy ONT
reads; a path distinguished by a single ONT base can also be a read
error rather than a biological RNA branch.
The [TSV](sv_transcript_evidence.tsv) retains original local witness
read IDs, raw alignment MAPQs, quality availability, intervals, and
assembly limitations. The collector allows missing base qualities;
no base-quality minimum is imposed on these counts.
Full-path exact support requires one read to match the entire candidate,
which is a much stricter observation than assembly from overlapping reads.
MAPQ 20 is a sensitivity subset, not a gate for the primary SV search.

| Sample / event | Path | RNA bases | Direct join (all / MAPQ20) | 12-flank exact | Selected local exact / window | Full exact | Limits |
|---|---|---:|---:|---:|---:|---:|---|
| HCC1395 chr13_52688329_chr13_60966095_+- | `ca1f6566` | 2430 | 4 / 3 | 1 (one retained) | 1 / 12 flank (one retained) | 1 | reference_models_unavailable |
| HCC1395 chr21_42975947_chr21_43035721_-+ | `be71afcf` | 2120 | 2 / 2 | 1 (unique) | 1 / 12 flank (unique) | 1 | reference_models_unavailable |
| HCC1395 chr21_42975947_chr21_43035721_-+ | `476eb4d7` | 2225 | 2 / 2 | 1 (unique) | 1 / 12 flank (unique) | 1 | reference_models_unavailable |
| HCC1395 chr9_83969153_chr9_84015113_-+ | `351071fd` | 1953 | 7 / 7 | 2 (shared) | 1 / 24 flank (unique) | 1 | extension_segment_limit, path_limit, reference_models_unavailable |
| HCC1395 chr9_83969153_chr9_84015113_-+ | `7c9ff8ab` | 1796 | 7 / 7 | 2 (shared) | 1 / 24 flank (unique) | 1 | extension_segment_limit, path_limit, reference_models_unavailable |
| HCC1395 chr9_83969153_chr9_84015113_-+ | `386f7562` | 2198 | 7 / 7 | 3 (shared) | 3 / 24 flank (unique) | 0 | extension_segment_limit, path_limit, reference_models_unavailable |
| COLO829 SV3_DEL | `88ab60e1` | 1213 | 3 / 0 | 3 (one retained) | 3 / 12 flank (one retained) | 0 | extension_segment_limit, reference_models_unavailable |

## Observed join versus unrearranged reference

The RNA column is `20 observed upstream bases | unplaced join bases |
20 observed downstream bases`. Reference columns are oriented genomic
controls at the observed left and right anchors: donor continuation
and acceptor predecessor, respectively. These are **genomic arms**,
not asserted mature reference transcript sequences. The full assembled
candidate sequences are in [`sv_transcript_evidence.fasta`](sv_transcript_evidence.fasta).

| Event / path | RNA join | Reference donor arm | Reference acceptor arm |
|---|---|---|---|
| chr13_52688329_chr13_60966095_+- `ca1f6566` | `CTTCTTCTTTAAGTCTTTTA|TCTT|CCCAGAGTATTAAAGCATTC` | `ACCTTTTTCTTAAGTTTTTA|ATGAGAACTTTGAGAAATTA` | `GTCACAATAGTCCAGCCTTC|CCCAGAGTATAAAAGCAACC` |
| chr21_42975947_chr21_43035721_-+ `be71afcf` | `TGGTCTCTACTAAAAATAAA|AAA|TGTGAACCTTGCCAGTATAG` | `TGGTCTCTACTAAAAATAAA|AAAATTAGTCGGGCATGGTG` | `AAGGAAGAATCGAGGTGAAA|TGTGAACCTCTGCCAGTATA` |
| chr21_42975947_chr21_43035721_-+ `476eb4d7` | `ACCTGGTCTACTAAAATAAA|AAA|TGTGAACCTTGCCAGTATAG` | `TGGTCTCTACTAAAAATAAA|AAAATTAGTCGGGCATGGTG` | `AAGGAAGAATCGAGGTGAAA|TGTGAACCTCTGCCAGTATA` |
| chr9_83969153_chr9_84015113_-+ `351071fd` | `GGCGTGCAAAAACTAACATT|GAAA|CAGTGAAAAAGAGGATTGCT` | `GGCGTGCAAAAACTAACATT|GAAATTTTGAAACAGCAGCA` | `TAGGAGTGGGACGTGGGAAA|CAGTGAAAAAGAGGATTGCT` |
| chr9_83969153_chr9_84015113_-+ `7c9ff8ab` | `GGCGTGCAAAAACTAACATT|GAAA|CAGTGAAAAAGAGGATTGCT` | `GGCGTGCAAAAACTAACATT|GAAATTTTGAAACAGCAGCA` | `TAGGAGTGGGACGTGGGAAA|CAGTGAAAAAGAGGATTGCT` |
| chr9_83969153_chr9_84015113_-+ `386f7562` | `GGGCGTGCAAAACTAACATT|GAAA|CAGTGAAAAAGAGGATTGCT` | `GGCGTGCAAAAACTAACATT|GAAATTTTGAAACAGCAGCA` | `TAGGAGTGGGACGTGGGAAA|CAGTGAAAAAGAGGATTGCT` |
| COLO829_SOMATIC_SV3_DEL `88ab60e1` | `AATTTGCCCTGTTGTCAGGC|CTATAGTAGAACTTCTTTT|GGAGTTGCTTCTTATGGATG` | `AATTTGCCCTGTTGTCAGGC|TGGAGTGCAGGGAAGTGATC` | `GATCAAATGTGAATGTTTGA|GGAGTTGCTTCTTATGGATG` |

The GRCh38 controls are pinned in
[`sv_reference_flanks.json`](../../tests/data/hcc1395_sv_junction/sv_reference_flanks.json)
from the [NCBI no-alt analysis set](https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.15_GRCh38/seqs_for_alignment_pipelines.ucsc_ids/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gz), uncompressed FASTA SHA-256 `9cce8b926416dd96b152deea85188495b75f7ac8d634cc723a017067be8702b7`.

## Remaining SV screen signals

The HCC1395 whole-run screen has six nominated SV records with
split-alignment candidates. The three outside this table are
chr1 records associated with POU2F1–NME7 or nearby complex
breaks; all three cite the **same single ONT read**
`SRR23456709.320596`. Its selected fixture contains two actual
alignment records and SA tags pointing to another placement;
an SA hint alone cannot reconstruct that missing piece. These
are one-read hypotheses, not three independent supported RNA
transcripts. The remaining COLO829 SV screen records have mate
context or breakend coverage but no independently observed
junction from which this method can assemble a path. The
[HCC1395](HCC1395_ATLAS.md) and [COLO829](COLO829_ATLAS.md)
source atlases enumerate the larger screened
call sets; this table covers all directly replayed joins with
at least two original fragments in the pinned selected-site
fixtures.

## Stronger fusion path from a separate HCC1395 source

The EIF3K–CYP39A1 fusion-enriched INTEGRATE RNA source supports one
531-aa predicted proteoform with 149 distinct Illumina junction
fragments, 19 tiling Illumina reads covering the 1,596-nt coding path,
and four independent ONT reads with exact junction sequence on
collinear start-to-stop alignments. Matched tumor DNA has 18 reciprocal
pairs between flanking introns versus zero normal pairs. These are
separate data and cannot be added to the GRCh38 selected-site counts.
See the [fusion report](eif3k_cyp39a1/report.md) and
[protein FASTA](eif3k_cyp39a1/candidate_proteins.fasta).

## What SV3 still lacks

SV3 has three original split-junction fragments at all MAPQs and
three exact 12-base-flank local witnesses, but zero MAPQ20 join
fragments and zero whole-path witnesses. Six opposite-mate pairs
are context, not split-junction support. The 1,213-base candidate
reaches `extension_segment_limit`, with no observed mature splice
connection to a CNIH3 exon. Its locus is intronic under Ensembl 87
and lacks a bundled ORFdata coding block. Two ATG-to-stop translations
(37 and 162 aa) share the same three nucleotide join witnesses;
neither has an established initiation site, reading frame, or full
ORF witness. A separate 68,553-representative polished isoform
check found no exact SV3 27-mer. Thus the RNA join is observed,
while a mature transcript, somatic causation and proteoform remain
unresolved. See the [SV3 ORF assessment](colo829_sv3_candidate_orfs.assessment.md).

The original BAMs are selected-site fixtures, not all HCC1395 or
COLO829 RNA reads. The HCC and COLO DNA/RNA sources also come from
different cell cultures or passages. The [HCC1395](HCC1395_ATLAS.md)
and [COLO829](COLO829_ATLAS.md) atlases give the screened source
universe and events without direct reconstructed joins.

Regenerate the table from pinned inputs with
`python -m reports.hcc1395_colo829.build_sv_transcript_evidence`.
The exact-window policy and its known ONT error sensitivity are
specified in `SV_TRANSCRIPT_EVIDENCE_SPEC.md`.

The [versioned exon path audit](sv_exon_path_evidence.md) is generated
by the same command and keeps exon connection and frame status
separate from local join evidence.
