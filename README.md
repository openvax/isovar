[![DOI](https://zenodo.org/badge/18834/hammerlab/isovar.svg)](https://zenodo.org/badge/latestdoi/18834/hammerlab/isovar) [![Build Status](https://travis-ci.org/hammerlab/isovar.svg?branch=master)](https://travis-ci.org/hammerlab/isovar)

# isovar
Abundance quantification of distinct transcript sequences containing somatic variants from cancer RNAseq

## Example

```sh
$ isovar-protein-sequences.py  \
    --vcf somatic-variants.vcf  \
    --bam rnaseq.bam \
    --genome hg19 \
    --min-reads 2 \
    --protein-sequence-length 30 \
    --output isovar-results.csv

  chr       pos ref alt                      amino_acids  \
0  22  46931060   A   C   FGVEAVDHGWPSMSSGSSWRASRGPPPPPR
1  22  46931062   G   A  CFGVEAVDHGWPPMSLAHGGPAVVHRLHPEA

   variant_aa_interval_start  variant_aa_interval_end ends_with_stop_codon  \
0                         16                       17                False
1                         16                       17                False

  frameshift  translations_count  supporting_variant_reads_count  \
0      False                   1                               1
1      False                   1                               1

   total_variant_reads  supporting_transcripts_count  total_transcripts  \
0                  130                             2                  2
1                  127                             2                  2

     gene
0  CELSR1
1  CELSR1
```

## Algorithm/Design

The one line explanation of isovar: `ProteinSequence = VariantSequence + ReferenceContext`.

A little more detail about the algorithm:
  1. Scan through an RNAseq BAM file and extract sequences overlapping a variant locus (represented by `ReadAtLocus`)
  2. Make sure that the read contains the variant allele and split its sequence into prefix/alt/suffix string parts (represented by `VariantRead`)
  3. Combine multiple `VariantRead` records into a `VariantSequence`
  4. Gather possible reading frames for distinct reference sequences around the variant locus (represented by `ReferenceContext`).
  5. Use the reading frame from a `ReferenceContext` to translate a `VariantSequence` into a protein fragment (represented by `Translation`).
  6. Multiple distinct variant sequences and reference contexts can generate the same translations, so we aggregate those equivalent `Translation` objects into a `ProteinSequence`.

Since we may not want to deal with *every* possible translation of *every* distinct sequence detected around a variant, `isovar` sorts the variant sequences by the number of supporting reads and the reference contexts in order of protein length and a configurable number of
translated protein fragments can be kept from this ordering.
