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

## Algorithm/Design

The one line explanation of isovar: `ProteinSequence = VariantSequence + ReferenceContext`.

A little more detail about the algorithm:
  1. Scan through an RNAseq BAM file and extract sequences overlapping a variant locus (represented by `ReadAtLocus`)
  2. Make sure that the read contains the variant allele and split its sequence into prefix/alt/suffix string parts (represented by `VariantRead`)
  3. Combine multiple `VariantRead` records into a `VariantSequence`
  4. Gather possible reading frames for distinct reference sequences around the variant locus (represented by `ReferenceContext`).
  5. Use the reading frame from a `ReferenceContext` to translate a `VariantSequence` into a protein fragment (represented by `ProteinSequence`).

Since we may not want to deal with *every* possible translation of *every* distinct sequence detected around a variant, `isovar` sorts the variant sequences by the number of supporting reads and the reference contexts in order of protein length and a configurable number of
translated protein fragments can be kept from this ordering.
