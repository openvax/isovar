# HCC1395 EIF3K-CYP39A1 fusion

[Open the focused PDF](HCC1395_EIF3K_CYP39A1_fusion.pdf).

**531-aa complete coding candidate**: 149 exact Illumina junction fragments; exact coding tiling includes 19 mapped reads. Separately, 9 ONT coding-span reads, 4 with an exact junction; 6/21 ONT candidates have an exact junction overall.

**IGV (GRCh37): `chr19:39123318-39123318 chr6:46607405-46607405`.**

EIF3K residues 1-166 are retained. Fusion residue 167 is a split codon (1 EIF3K base + 2 CYP39A1 bases) encoding the same alanine as CYP39A1 residue 105; the remaining sequence matches CYP39A1 residues 106-469. The joined context is novel, not the entire suffix.

D = retained EIF3K; F = fusion; A = retained CYP39A1. Dashes mean a contributor is absent.

```text
D MAMFEQMRANVGKLLKGIDRYNPENLATLERYVETQAKENAYDLEANLAVLKLYQFNPAF
F MAMFEQMRANVGKLLKGIDRYNPENLATLERYVETQAKENAYDLEANLAVLKLYQFNPAF
A ------------------------------------------------------------

D FQTTVTAQILLKALTNLPHTDFTLCKCMIDQAHQEERPIRQILYLGDLLETCHFQAFWQA
F FQTTVTAQILLKALTNLPHTDFTLCKCMIDQAHQEERPIRQILYLGDLLETCHFQAFWQA
A ------------------------------------------------------------

D LDENMDLLEGITGFEDSVRKFICHVVGITYQHIDRWLLAEMLGDLS--------------
F LDENMDLLEGITGFEDSVRKFICHVVGITYQHIDRWLLAEMLGDLSASIPKNVFLALHEK
A ----------------------------------------------ASIPKNVFLALHEK

D ------------------------------------------------------------
F LYIMLKGKMGTVNLHQFTGQLTEELHEQLENLGTHGTMDLNNLVRHLLYPVTVNMLFNKS
A LYIMLKGKMGTVNLHQFTGQLTEELHEQLENLGTHGTMDLNNLVRHLLYPVTVNMLFNKS

D ------------------------------------------------------------
F LFSTNKKKIKEFHQYFQVYDEDFEYGSQLPECLLRNWSKSKKWFLELFEKNIPDIKACKS
A LFSTNKKKIKEFHQYFQVYDEDFEYGSQLPECLLRNWSKSKKWFLELFEKNIPDIKACKS

D ------------------------------------------------------------
F AKDNSMTLLQATLDIVETETSKENSPNYGLLLLWASLSNAVPVAFWTLAYVLSHPDIHKA
A AKDNSMTLLQATLDIVETETSKENSPNYGLLLLWASLSNAVPVAFWTLAYVLSHPDIHKA

D ------------------------------------------------------------
F IMEGISSVFGKAGKDKIKVSEDDLENLLLIKWCVLETIRLKAPGVITRKVVKPVEILNYI
A IMEGISSVFGKAGKDKIKVSEDDLENLLLIKWCVLETIRLKAPGVITRKVVKPVEILNYI

D ------------------------------------------------------------
F IPSGDLLMLSPFWLHRNPKYFPEPELFKPERWKKANLEKHSFLDCFMAFGSGKFQCPARW
A IPSGDLLMLSPFWLHRNPKYFPEPELFKPERWKKANLEKHSFLDCFMAFGSGKFQCPARW

D ---------------------------------------------------
F FALLEVQMCIILILYKYDCSLLDPLPKQSYLHLVGVPQPEGQCRIEYKQRI
A FALLEVQMCIILILYKYDCSLLDPLPKQSYLHLVGVPQPEGQCRIEYKQRI

Terminal stop: *
```

The PDF includes both complete reference proteins, a per-codon exact-fragment coverage track, and an assembly track. All 1,596 coding bases are RNA-witnessed; no bases are filled from reference. The junction is highlighted separately from retained parental sequence. No additional germline or phased somatic amino-acid change is observed.

Sources: INTEGRATE public example RNA and tumor/normal DNA, independent SRR23456709 ONT RNA, archive GRCh37-lite transcripts with sequence-verified Ensembl 87 CDS, and separate GRCh38 ONT exon models. See the [complete source index](../SEQUENCING_SOURCES.md) and adjacent JSON for every input hash and per-codon exact-fragment support.
