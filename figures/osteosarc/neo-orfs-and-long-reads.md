# Sid non-SNV neo-ORFs and long-read evidence

Reviewed 2026-09-21 with Isovar 1.21.4 and Varcode 9.3.7. Varcode was upgraded in the shared environment; all 4,564 Isovar tests and lint passed. The counts below separate full-locus historical audits from selected regression fixtures. Processing products are not pooled.

## Main finding: TPST1–CRCP has a well-supported candidate upstream ORF

The previous phrase “no novel ORF” was too strong. The annotation-based pipeline reports no donor-CDS-anchored translation across this breakpoint, but an ATG scan of the observed RNA identifies this candidate:

`MPSRRRGGSSLIPGRMPILRFSVTTRYFSY*`

This is 30 amino acids plus a stop. Its complete 93-nt ATG-to-stop sequence is observed in the original stored reads, not filled from reference:

| Product | Exact full-sequence support | Quality and counting limits |
|---|---:|---|
| T1 ONT tagged, selected junction fixture | 96 distinct RG/QNAMEs; 61 RG/CB/UB labels | Q10 or better at every candidate base |
| T2 ONT tagged, selected junction fixture | 132 distinct RG/QNAMEs; 122 RG/CB/UB labels | Q10 or better at every candidate base |
| T1 PacBio, original breakpoint-path fixture | 15 distinct RG/QNAMEs and RG/CB/XM labels | QUAL unavailable; all retain the actual supplementary partner |

These are fixture counts, not whole-library abundance estimates or proof of independent molecules. T1 ONT and PacBio share the single-cell cDNA source and are not independent biological replication. Of the 15 PacBio witnesses, 13 have `rc=1`, two `rc=0`; this is a cell-calling annotation, not a reason to discard RNA sequence or assign malignancy.

The candidate ATG lies at GRCh38 chr7:66205483 (1-based), in the first exon/5′ UTR of TPST1 ENST00000304842. It is at cDNA offset 284, 141 nt before the annotated CDS start at offset 425. The fusion leaves TPST1 before that annotated start. Thus there is no annotated TPST1 coding frame to transfer across the junction. The acceptor's annotated frame cannot establish initiation upstream on the donor.

The annotated TPST1 start is in exon 2 at chr7:66240426 (1-based). The fusion joins exon 1 to CRCP instead of continuing into that exon. Using the normal start would insert sequence absent from the observed RNA.

The start context is `GCCAGG-ATG-CCGTCC`: A at −3 but C at +4. We enumerated ATGs; we did not select the first ATG, optimize Kozak scores, or assume that the longest ORF is translated. In the reference TPST1 transcript, the same unannotated ATG yields a different potential **39-aa** upstream ORF (`MPSRRRGGSSLIPDVGYLSEVDCPWPEHFPKIILSKISV`). The earlier prose incorrectly counted it as 37 aa; the reconstructed sequence was correct. The fusion changes its continuation after the shared 13-aa prefix. This is a **candidate altered upstream ORF**, not evidence that initiation or protein expression occurs in Sid. Genome-wide peptide novelty has not been checked.

[Kowar et al. (2025)](https://doi.org/10.1038/s41467-025-63405-2) independently
include this exact 39-aa sequence in their immunopeptidomics search database
(Supplementary Data 12, row 2227; its trailing `Z` is not an amino acid in this
ORF). Search-database inclusion is not an observed MS peptide. Their Ribo-seq
tables report longer TPST1 upstream ORFs in the same frame: an ACG start in
nocodazole- and Taxol-treated U-2 OS cells (Supplementary Data 2 and 4), and a
CTG start in nocodazole-treated MDA-MB-231 cells (Data 9). These calls support
translation of the upstream region under those conditions but do not establish
initiation at our internal ATG or translation of Sid's fusion.

We searched the study's three processed immunopeptidomics tables in
[PRIDE PXD057839](https://www.ebi.ac.uk/pride/archive/projects/PXD057839): U-2 OS
(15,880 precursor rows), SUM-159PT (30,314) and MDA-MB-231 (20,294). None reported
TPST1/ENST00000304842 or a peptide of at least seven residues contained in the
39-aa wild-type sequence, including a search treating isoleucine/leucine as
indistinguishable. These are precursor rows, not distinct peptides or biological
replicates. Thus this study provides regional Ribo-seq support, but no detected
MS peptide from this uORF in its processed results. That is not proof of absence
of translation in another tissue, condition or study.

Alternative PacBio paths contain 61- and 72-aa candidates starting at the same ATG; each has only one exact full ATG-to-stop witness. The 20-template breakpoint support must not be assigned to every complete candidate. Their differences remain possible splice/error/assembly alternatives.

Long reads demonstrate full candidate sequence linkage, but the 93-nt ORF also fits inside the previously retained 120-nt junction window. We have not established that short reads cannot resolve it or that long reads uniquely determine initiation.

## Other non-SNV candidates

| Candidate | Evidence | Reading-frame / long-read conclusion |
|---|---|---|
| **TECPR1**, 1-nt deletion, chr7:98241127 TG>T | Current selected T2 Illumina fixture: 17 alt templates; 15 contribute to a 49-aa assembly with 25 aa in the shifted frame. Without assembly: 37 aa, with the same 25-aa shifted tail. | Strong frameshift example exceeding 10 contributing templates. The extra assembled context is upstream reference context. No single read spans the full 49-aa assembly; no stop is observed. This does not demonstrate a long-read advantage. |
| **PIP5K1A**, 1-nt deletion | Current T1 ONT: 4 alt templates, 2 exact full-window witnesses; 46-aa reconstruction, including 21 shifted aa and an observed stop. Matched-timepoint Illumina: 1 alt template, no protein. | Longer linked context exists in ONT, but support is below >10 and unequal alternate depth prevents attributing the result solely to read length. Frame transferred from annotated CDS. |
| **GLIS3**, 8-nt deletion | Historical full-locus products reach 19 supporting read objects but only 10 template names. A T1 short-read reconstruction has 25 shifted aa. | Does not exceed 10 distinct templates; no ONT alternate support in the bounded comparison. |
| **KTN1**, 2-nt insertion | T2 ONT: 10 alt templates, 9 contribute to a 29-aa protein; short reads: 3 templates, 18 aa. Both recover 5 shifted aa. | ONT adds upstream context, not a longer novel tail; below >10 template support. |
| **FOXO3–STRADA/CCDC47** | Current tagged ONT breakpoint path: 12 fragments / 8 molecule labels. | Frame unresolved; no justified protein from the present coding-anchor method. Candidate alternative initiation needs a separate audit. |
| **ATP5MG–KMT2A** | Short-read fixture: two exact witnesses for an 18-aa conditional protein including an observed annotated ATP5MG start and stop. PacBio: five direct observations of that join. | Below threshold; the same join is compatible with ordinary read-through splicing. Other longer PacBio protein hypotheses are not the same nominated junction and have weak support. |
| **GTF3C5**, **RNF213**, **H1-2** deletions | Numerous alternate reads; GTF3C5 has 117/128 T1/T2 ONT alt templates; RNF213 has 42 T1 short-read alt templates. | In-frame deletions, not long frameshift neo-ORFs. Do not equate total protein-window length with new sequence length. |

CD109 is still a useful long-read phasing example (21 templates carrying two SNVs), and ZNF436 a context-length example, but neither is a non-SNV neo-ORF example.

## PacBio quality and provenance

All 522 bundled PacBio records were inspected: 515 lack QUAL. Their header records `READTYPE=TRANSCRIPT`, Iso-Seq `groupdedup` 4.0.0 and pbmm2 `--preset ISOSEQ`. Missing QUAL is a property of this processed product, not of every aligned PacBio read.

Observed tags: `CB`, `NM`, `RG`, `SA`, `XA`, `XM`, `ic`, `im`, `is`, `it`, `mg`, `rc`, `rm`, `zm`. No `rq`, `np`, `ec`, or `sn` tags are present. `ic/is` describe inputs to the transcript consensus; they are not CCS polymerase pass counts and are not extra independent molecules. `im` preserves input names. `CB/XM` preserve barcode/UMI provenance. `NM`, `mg`, and MAPQ describe alignment evidence, not a calibrated probability for each candidate base. Unknown tags are retained without guessing their semantics.

`rm=1` specifically records pbmm2's trimming of overlapping query matches. `SA` lists supplementary alignments, `RG` the read group, `XA` the Iso-Seq tag-name order, `it` the clipped barcode/UMI list, and `zm` the ZMW identifier. Tags must be interpreted with their producer: Iso-Seq's `XA` is not another aligner's alternative-alignment tag.

Keep missing base qualities optional (already the default). Expose the available sequence, mapping, consensus and provenance evidence separately; do not synthesize Phred scores. For the 15 exact PacBio candidate witnesses, 14 have `ic=is=1`, one has `ic=is=2`, and they retain 16 input CCS names. Exact ONT recurrence supplies additional sequence evidence with measured qualities.

## Sources and reproducibility

The production exploratory ORF output and all supporting read tags can now be
reproduced offline from the package's existing minimal read bundle:

```sh
python -m examples.osteosarc_sv_orfs /tmp/tpst1-orfs.json
```

`result.paths[].exploratory_orfs` contains the candidates and full-span witness
IDs; `result.record_evidence` and `result.original_records` expose their native
tags and SAM records. The original annotation-driven translations remain
separate. This command reproduces the PacBio 30/61/72-aa comparison; the ONT
counts above refer to the separate selected corpus fixtures.

- Primary data and annotation: `tests/data/fusions/corpus/TPST1--CRCP-{T1,T2}.input.json.gz`, `tests/data/fusions/long-read/TPST1--CRCP.PacBio-T1.sam.gz`, and their checksum manifests.
- Current reanalysis scripts/results: `/private/tmp/isovar-current-neoorf-audit.py`, `/private/tmp/isovar-current-neoorf-audit.json`, `/private/tmp/inspect-tpst1-orf.py`, `/private/tmp/tpst1-current-full.json`, `/private/tmp/tpst1-candidate-support.py`, `/private/tmp/tpst1-candidate-support.json`, `/private/tmp/isovar-pacbio-tag-inventory.json`.
- Historical counts: `tests/data/osteosarc/expansion/audit/matrix.json.gz`; `tests/data/osteosarc/figure_comparisons/{README,RNA_FOOTPRINTS,EXTENDED_RNA}.md`.
- [TECPR1 source allele](https://osteosarc.com/variant/TECPR1-chr7-98241127/); [fusion catalogue](https://osteosarc.com/fusions/).
- [PacBio BAM specification](https://pacbiofileformats.readthedocs.io/en/13.1/BAM.html), [Iso-Seq tag definitions](https://isoseq.how/isoseq-tags.html), [pbmm2 documentation](https://github.com/PacificBiosciences/pbmm2).
- [Kozak's initiation-context mutagenesis](https://pubmed.ncbi.nlm.nih.gov/3943125/); [NCBI ORFfinder](https://www.ncbi.nlm.nih.gov/orffinder/).

The follow-up implementation for [#326](https://github.com/openvax/isovar/issues/326) and [#327](https://github.com/openvax/isovar/issues/327) adds explicit exploratory ORF candidates, annotation/start-context labels, full-interval support and structured native read tags. Annotated-CDS translations remain distinct. Long/short comparisons still must not combine processing copies or infer absent expression from absent support.
