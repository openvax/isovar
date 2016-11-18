[![DOI](https://zenodo.org/badge/18834/hammerlab/isovar.svg)](https://zenodo.org/badge/latestdoi/18834/hammerlab/isovar) [![Build Status](https://travis-ci.org/hammerlab/isovar.svg?branch=master)](https://travis-ci.org/hammerlab/isovar) [![Coverage Status](https://coveralls.io/repos/github/hammerlab/isovar/badge.svg?branch=master)](https://coveralls.io/github/hammerlab/isovar?branch=master)

# isovar
Isovar assembles protein subsequences around mutations from cancer RNA-Seq data. Since Isovar uses sequenced reads to determine a mutant coding sequence it is able to correctly phase somatic variants with adjacent germline variants, as well as sometimes recovering alternatively spliced isoforms. 

## Example

```sh
$ isovar-protein-sequences.py  \
    --vcf somatic-variants.vcf  \
    --bam rnaseq.bam \
    --min-reads 2 \
    --protein-sequence-length 30 \
    --output isovar-results.csv

  chr       pos ref alt                      amino_acids  \
0  22  46931060   A   C   FGVEAVDHGWPSMSSGSSWRASRGPPPPPR

   variant_aa_interval_start  variant_aa_interval_end ends_with_stop_codon  \
0                         16                       17                False

  frameshift  translations_count  supporting_variant_reads_count  \
0      False                   1                               1

   total_variant_reads  supporting_transcripts_count  total_transcripts     gene
0                  130                             2                  2   CELSR1
```

## Algorithm/Design

The one line explanation of Isovar: `ProteinSequence = VariantSequence + ReferenceContext`.

A little more detail about the algorithm:
  1. Scan through an RNAseq BAM file and extract sequences overlapping a variant locus (represented by `LocusRead`)
  2. Make sure that the read contains the variant allele and split its sequence into prefix/alt/suffix string parts (represented by `AlleleRead`)
  3. Assemble overlapping `AlleleRead`s (which agree with the variant allele) into a `VariantSequence`
  4. Gather possible reading frames for distinct reference sequences around the variant locus (represented by `ReferenceContext`).
  5. Use the reading frame from a `ReferenceContext` to translate a `VariantSequence` into a protein fragment (represented by `Translation`).
  6. Multiple distinct variant sequences and reference contexts can generate the same translations, so we aggregate those equivalent `Translation` objects into a `ProteinSequence`.

Since we may not want to deal with *every* possible translation of *every* distinct sequence detected around a variant, Isovar sorts the variant sequences by the number of supporting reads and the reference contexts in order of protein length and a configurable number of translated protein fragments can be kept from this ordering.

## Sequencing Recommendations

Isovar works best with high quality / high coverage mRNA sequence data. This means that you will get best results from >=100M paired-end reads sequenced on an Illumina HiSeq from a library enriched with poly-A capture. The length of 
