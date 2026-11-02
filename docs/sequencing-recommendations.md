# Sequencing recommendations

This page is for anyone planning or choosing RNA-seq data for Isovar: what
kind of library works best, and how read length and platform limit the protein
it can reconstruct.

Isovar works best with high-quality, high-coverage poly-A-selected mRNA
sequencing, for example >100M paired-end reads on a current Illumina short-read
platform. The depth needed depends on RNA degradation and tumor purity.

With short reads, read length bounds the recoverable protein. Assembly only uses
reads overlapping the variant, so 100 bp reads give at most 199 bp of sequence
around a somatic SNV, about 66 amino acids. Without assembly, one 100 bp read
determines at most 33.

Overlap assembly requires exact sequence matches, which suits short reads with
low error rates. Long reads (PacBio, Oxford Nanopore) often span the whole
context without assembly, but noisy reads may not join by exact overlap; SV
reconstruction (`isovar sv-rna`) uses noise-tolerant extension. Coverage
trimming assumes that read coverage falls off away from the variant, which reads
spanning splice junctions can violate.

See [read eligibility across platforms](read-processing-audit.md) for which
alignments Isovar uses from Illumina, ONT and PacBio data.
