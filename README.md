<a href="https://travis-ci.org/openvax/isovar">
    <img src="https://travis-ci.org/openvax/isovar.svg?branch=master" alt="Build Status" />
</a>
<a href="https://coveralls.io/github/openvax/isovar?branch=master">
    <img src="https://coveralls.io/repos/openvax/isovar/badge.svg?branch=master&service=github" alt="Coverage Status" />
</a>
<a href="https://pypi.python.org/pypi/isovar/">
    <img src="https://img.shields.io/pypi/v/isovar.svg?maxAge=1000" alt="PyPI" />
</a>

# isovar
Isovar assembles protein subsequences around mutations from cancer RNA-Seq data. Since Isovar uses sequenced reads to determine a mutant coding sequence it is able to correctly phase somatic variants with adjacent germline variants, as well as sometimes recovering alternatively spliced isoforms.

## Example

```sh
$ isovar-protein-sequences  \
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

Isovar works best with high quality / high coverage mRNA sequence data. This means that you will get best results from >100M paired-end reads sequenced on an Illumina HiSeq from a library enriched with poly-A capture. The number of reads varies depending on degree of RNA degradation and tumor purity. The read length will determine the longest protein sequence you can recover, since Isovar's cDNA assembly only considers reads that overlap a variant. With 100bp reads you will be able to assemble at most 199bp of sequence around a somatic single nucleotide variant, and consequently only be to determine 66 amino acids from the protein sequence. If you disable the cDNA assembly algorithm then a 100bp read will only be able to determine 33 amino acids.

## Commandline Tools

* `isovar-protein-sequences`
* `isovar-allele-counts`
* `isovar-allele-reads`
* `isovar-translations`
* `isovar-reference-contexts`
* `isovar-variant-reads`
* `isovar-variant-sequences`

