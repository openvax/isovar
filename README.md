[![DOI](https://zenodo.org/badge/18834/hammerlab/isovar.svg)](https://zenodo.org/badge/latestdoi/18834/hammerlab/isovar)

# isovar
Abundance quantification of distinct transcript sequences containing somatic variants from cancer RNAseq

## Example

```sh
$ isovar-translate-variants.py \
    --vcf /Users/iskander/code/varlens/test/data/CELSR1/vcfs/vcf_1.vcf  \
    --genome hg19 \
    --bam /Users/iskander/code/varlens/test/data/CELSR1/bams/bam_2.bam  \
    --protein-fragment-length 10 \
    --min-reads 20 \
    --output isovar-results.csv

Variant(contig=22, start=46931062, ref=G, alt=A, reference_name=GRCh37)
  chr  base1_start_pos  base1_end_pos ref alt variant_protein_sequence  \
0  22         46931060       46931060   A   C               PPMSSATSVS
1  22         46931062       46931062   G   A              WPPMSFSTSVS

   variant_protein_sequence_length         reference_transcript_ids  \
0                               10  ENST00000395964;ENST00000262738
1                               11  ENST00000395964;ENST00000262738

  reference_transcript_names                       cdna_sequences  \
0      CELSR1-201;CELSR1-001  GCCCCCCATGAGCTCC_G_CCACCAGCGTGTCCAT
1      CELSR1-201;CELSR1-001  TGGCCCCCCATGAGCT_T_CTCCACCAGCGTGTCC

   cdna_sequence_length  number_supporting_reads
0                    33                       48
1                    33                       49
```

