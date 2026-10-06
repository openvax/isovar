# Sid SV-associated and frameshift RNA candidates

The focused October 5 search identifies **TPST1–CRCP, SPAG1, TECPR1,
PIP5K1A and FOXO3–STRADA/CCDC47** as candidates worth following up.
SPAG1 is an additional non-vaccine frameshift lead. The evidence below is
original RNA sequence covering each stated interval. Translation, tumor
specificity and antigen presentation have not been demonstrated.

The [complete DNA → RNA → ORF report](sid-neoorf-event-report.md) adds original
DNA calls/alleles and support, transcript/event anatomy, every sequence and
per-product support, including unresolved and zero-candidate entries.

## Main candidates

Fragment counts refer to complete coverage of the stated sequence window,
deduplicated by read group and query name. Q20 means every base in that interval
has a recorded quality of at least 20. A complete coding window need not cover
the complete gene ORF or its initiation site. Products and timepoints are kept
separate; shorter and longer windows often share witnesses.

| Candidate | RNA evidence | Interpretation and remaining uncertainty |
| --- | --- | --- |
| **TPST1–CRCP**, altered upstream ORF | Exact 30-aa sequence plus stop: **172 ONT T2 fragments**, 139 barcode labels; **15 PacBio T1 fragments**, 15 labels. | Best supported altered-uORF candidate in the selected joins. Its unannotated ATG is upstream of the TPST1 CDS; the fusion changes the continuation after a shared 13-aa prefix. Initiation and RNA polarity remain unresolved. PacBio qualities are unavailable; SV interval quality was not assessed. |
| **SPAG1 chr8:100213209**, frameshift | 10x T2: **16 fragments, 15 Q20** cover a window with 10 shifted aa; **2 fragments, both Q20** cover a longer **25-aa shifted tail**. There are 20 alternate-allele fragments at the locus. | Strong new non-vaccine lead with an annotated coding frame. The longer tail is supported by fewer reads; no stop or initiation site is established by these local windows. |
| **TECPR1 chr7:98241127**, frameshift | Bulk SARC0277: **15 fragments, 14 Q20** cover an 11-aa shifted tail; **2 fragments, both Q20** cover a **39-aa shifted tail**. There are 17 alternate-allele fragments. | Longer novel coding continuation is present in native reads, although its full-window support is two. No stop is observed in these windows. This locus was already included in the vaccine input panel. |
| **PIP5K1A chr1:151242178**, frameshift | ONT T2: **2 complete fragments, 1 Q20**, cover **21 shifted aa plus a stop**. 10x T2: **2 fragments, both Q20**, cover 17 shifted aa. Bulk has shorter windows with lower quality support. | Annotated coding frame and a linked stop in ONT. The ONT and 10x windows differ and their counts must not be pooled. The ONT locus has seven alternate-allele fragments, exceeding the two complete-window witnesses. |
| **FOXO3–STRADA/CCDC47**, exploratory ORFs | A representative **34-aa sequence plus stop** has **5 ONT T2 fragments**, five labels, and **3 PacBio T1 fragments**, three labels. Nested alternative starts also yield 36- and 67-aa hypotheses. | Original RNA supports sequence across the nominated breakpoint, but multiple frames/starts remain unresolved. Support is reverse-complement-only for these hypotheses and RNA polarity is unresolved. This does not establish which ORF is translated. |

The TPST1 sequence is `MPSRRRGGSSLIPGRMPILRFSVTTRYFSY*`.
The representative FOXO3 sequence is
`MEKLCSLYYGCVPLLLRYISITEYPYNNGILLSN*`.
SPAG1's longer window is `GAPQRARPRRPARTSGAHGGPLRRRRRAA`, with the
shifted interval beginning at zero-based amino-acid offset 4. TECPR1's longer
window is `PHDLLWAHSGRDRPWSGKESTGAIPKEVPGPSWSLLDLKTGSCTS`, with the shifted
interval beginning at offset 6. All of these sequences contain 8-mers absent
from the pinned annotated proteome after treating I and L as equivalent.

The earlier [TPST1 annotation analysis](../figures/osteosarc/neo-orfs-and-long-reads.md)
locates its candidate ATG in exon 1/5′ UTR, before the annotated CDS, and compares
it with the different potential 39-aa reference uORF. Regional Ribo-seq evidence
in [Kowar et al.](https://doi.org/10.1038/s41467-025-63405-2) does not demonstrate
initiation at this exact internal ATG or translation of Sid's fusion sequence.

## Other results

KTN1 chr14:55627965 has an annotated frameshift window with **11 ONT T2
fragments, eight Q20**, covering five shifted aa and a stop. It is a short
altered tail rather than a long new ORF. RNF213, H1_2 and GTF3C5 contribute
in-frame coding controls and are not labeled frameshift neoORFs.

PARD3B–CDKN2B-AS1 has a 10-aa exploratory hypothesis supported by 11 ONT T2
fragments and five PacBio T1 fragments. The two nominated geometries reuse
witnesses and must not be added. OTUD7A–FMN1 has a 21-aa hypothesis supported
by nine ONT fragments across seven labels. GABBR1–SLC29A1 has a 67-aa PacBio
hypothesis with nine fragments across eight labels, but its join is only
event-compatible. These all retain unresolved initiation/frame/polarity flags.

The dense ::IGKC screen has a 27-aa splice-ambiguous hypothesis with 279
fragments across 12 labels. Immune-receptor sequence and splice ambiguity make
it unsuitable to call a tumor-specific neoORF from these counts. ATP5MG–KMT2A
produces 29 splice-ambiguous hypotheses with **zero complete-ORF witnesses**;
its RNA join is compatible with ordinary read-through splicing. This is not a
negative claim about expression or all possible ATP5MG ORFs.

## Search scope and verification

The screen covers **nine selected SV geometries in 60 oriented views**: seven
geometries across four products, plus the earlier CNN2 and ::IGKC dense 10x
controls. It also screens **31 literal indel/splice nominations across four
products (124 allele/product screens)**. All 124 completed without screen
errors. MUC3A's nonliteral duplication is separately unassessable in four
products. Splice-site nominations are inspected as focal alleles here; this
does not discover arbitrary altered splicing, retained introns or exonization.

The [candidate ledger](data/sid-neoorf/candidate-screen.json.gz) contains
**1,302 source-specific hypotheses**: 412 SV RNA ORFs, 116 frameshift coding
windows and 774 in-frame windows. Overlapping windows, alternative starts,
frames and geometry aliases mean these are not 1,302 independent neoORFs.
It retains zero-result allele screens, acquisition bounds, RNA evidence,
transcript IDs, sequence identities and every uncertainty flag.

The comparison scans **245,535 Ensembl 115 protein entries** (116,853,147 aa),
with exact and I/L-collapsed 8-mer comparisons and separate full-sequence hits.
The protein FASTA SHA-256 is
`7991ab37d09986eefdb5c77032b15370ab0e99918b78f235f3230c0302851c96`.
Absence from these annotated proteins does not establish absence from normal
noncoding ORFs, matched normal tissue or the patient's other alleles. No HLA
binding or presentation screen is included. Consequences and reading frames
are transcript-specific, as described in the
[Ensembl consequence definitions](https://www.ensembl.org/info/genome/variation/prediction/predicted_data.html).

The [native coding-window checks](data/sid-neoorf/coding-window-checks.json.gz)
independently translate original SAM sequence in both orientations and confirm
the nucleotide hash, fragment identities and interval quality for nine selected
frameshift windows. Their original SAM record hashes are preserved. The
addition ledger independently verifies the FOXO3 representative against 13
original BAM records; prior ledgers verify the TPST1 and dense representatives.
Missing SAM qualities remain unavailable, consistent with the
[SAM specification](https://samtools.github.io/hts-specs/SAMv1.pdf).

## Reproduce

Install `isovar[data]` and the pinned Ensembl 115 annotation. Starting with the
[focused archive](sid-priority-read-subsets.md),
[additional regional inputs](sid-neoorf-read-subsets.md) and corrected dense
audit, run:

```sh
python -m examples.screen_sid_non_snv_orfs PRIOR_AUDIT VARIANT_RESULTS
python -m examples.screen_sid_non_snv_orfs ADDITIONS VARIANT_RESULTS
python -m examples.screen_sid_orf_candidates SCREEN.json.gz \
  --audit PRIOR_AUDIT --audit ADDITIONS --audit DENSE_AUDIT \
  --variants VARIANT_RESULTS --proteome ENSEMBL_115_PEP_FASTA
python -m examples.check_sid_coding_windows SCREEN.json.gz CHECKS.json.gz \
  --audit PRIOR_AUDIT --audit ADDITIONS
```

Checkpoints pin code, engine/settings, BAM, selection, Ensembl source files
and consumed annotation indexes. The release ledger preserves the successful
reference-pinned pass; earlier incomplete/annotation-alias attempts remain in
the local working cache and are excluded. Use fresh output paths after changes.
Relative archive pins verify copied inputs, while resumable receipts still use
original paths; portable replay remains tracked in
[#440](https://github.com/openvax/isovar/issues/440).
The complete 1,493-geometry catalogue and general neoORF discovery remain
outside this focused search, not negative results.
