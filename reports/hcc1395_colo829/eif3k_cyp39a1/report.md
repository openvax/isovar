# HCC1395 EIF3K–CYP39A1 coding fusion candidate

Open the [full protein alignment and per-codon RNA support PDF](../focused_proteins/HCC1395_EIF3K_CYP39A1_fusion.pdf). See the [complete sequencing and reference source index](../SEQUENCING_SOURCES.md).

- Original fusion-enriched Illumina records screened: 185
- Exact junction-spanning fragments: 149
- Original mapped Illumina reads in a coding-region overlap chain: 19
- Fully observed coding RNA path: 1596 bases
- Original ONT candidate reads with exact 20+20-base junction sequence: 6/21
- ONT reads spanning coding start to stop in collinear local alignments: 9/21
- Those with an exact junction sequence: 4 (1513–1539 aligned coding bases each)
- Tumor/normal DNA fragments linking the flanking introns: 18/0
- Isovar frame status: translated
- Candidate protein lengths: 531 amino acids
- Fusion junction: coding nucleotide 499, within amino-acid codon 167
- Example junction-spanning 9-mers: AEMLGDLSA, EMLGDLSAS, MLGDLSASI, LGDLSASIP, GDLSASIPK

Every base of the coding path is witnessed by an exact original read, and the selected reads form a chain with at least 20 bases of overlap. The exon path is reference guided; six distinct original ONT reads also contain the exact junction sequence. For this supplied coding path, Isovar finds one annotated frame with an observed start and stop. Alternative RNA isoforms remain possible. The tumor-specific cross-chromosome DNA pairs support a rearrangement between the flanking introns, without resolving its exact breakpoint. These data do not measure translated protein expression. See candidate_proteins.fasta for the sequence.

[Five-panel fusion figure](figures.pdf) shows the junction, annotated exon context, protein and example peptides.

The Illumina input is a fusion-enriched subset, so these counts are not whole-sample prevalence estimates. The ONT candidates came from a local two-chromosome alignment and are checked here by exact original query sequence at the join. The [INTEGRATE public example](https://sourceforge.net/p/integrate-fusion/wiki/example/) supplies the Illumina RNA and tumor/normal DNA sources; its [primary publication](https://pmc.ncbi.nlm.nih.gov/articles/PMC4691743/) describes the multiple EIF3K–CYP39A1 isoforms. The independent ONT run is [SRR23456709](https://www.ncbi.nlm.nih.gov/sra/SRR23456709).
