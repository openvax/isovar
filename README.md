<a href="https://travis-ci.org/openvax/isovar">
    <img src="https://travis-ci.org/openvax/isovar.svg?branch=master" alt="Build Status" />
</a>
<a href="https://coveralls.io/github/openvax/isovar?branch=master">
    <img src="https://coveralls.io/repos/openvax/isovar/badge.svg?branch=master&service=github" alt="Coverage Status" />
</a>
<a href="https://pypi.python.org/pypi/isovar/">
    <img src="https://img.shields.io/pypi/v/isovar.svg?maxAge=1000" alt="PyPI" />
</a>

# Isovar
Isovar determines mutant protein subsequences around mutations from cancer RNAseq data.

Isovar works by:

(1) collecting RNA reads which spanning the location of a variant,

(2) filtering the RNA reads to those which support the mutation,

(3) assembling mutant reads into longer coding sequences,  

(4) matching mutant coding sequences against reference annotated reading
frames, and

(5) translating coding sequences determined directly from RNA into mutant protein sequences.

The assembled coding sequences may incorporate proximal 
(germline and somatic) variants, along with any splicing alterations 
which occur due to modified splice signals.


## CLI Example

```sh
$ isovar  \
    --vcf somatic-variants.vcf  \
    --bam rnaseq.bam \
    --protein-sequence-length 30 \
    --output isovar-results.csv
```

## Python Example

```python

from isovar import run_isovar

isovar_results = run_isovar(
    variants="cancer-mutations.vcf",
    alignment_file="")
    
# this code traverses every variant and prints the number
# of RNA reads which support the alt allele for variants
# which had a successfully assembled/translated protein sequence
for isovar_result in isovar_results:
    if isovar_result.top_protein_sequence is not None:
        print(isovar_result.num_alt_fragments)
```

## Organization and Algorithmic Design

The inputs to Isovar are one or more somatic variant call (VCF) files, along with a BAM file 
containing aligned tumor RNA reads. 

* [LocusRead](https://github.com/openvax/isovar/blob/master/isovar/locus_read.py): Isovar examines each variant locus and extracts reads overlapping that locus, 
represented by `LocusRead`. The `LocusRead` representation allows filtering  based
on quality and alignment criteria (e.g. MAPQ > 0) which are thrown away in later stages
of Isovar. 

* [AlleleRead](https://github.com/openvax/isovar/blob/master/isovar/allele_read.py): Once `LocusRead` objects have been filtered, they are converted into a simplified 
representation called `AlleleRead`. Each `AlleleRead` contains only the cDNA sequences 
*before*, *at*, and *after* the variant locus. 

* [ReadEvidence](https://github.com/openvax/isovar/blob/master/isovar/read_evidence.py): 
The set of `AlleleRead` objects overlapping a mutation's location may support many different
distinct allele. The `ReadEvidence` type represents the grouping of these reads into
*ref*, *alt* and *other* `AlleleRead` sets, where *ref* reads agree with the reference
 sequence, *alt* reads agree with the given mutation, and *other* reads contain all
 non-ref/non-alt alleles. The *alt* reads will be used later to determine
a mutant coding sequence, but the *ref* and *other* groups are also kept in case they are
useful for filtering. 

* [VariantSequence](https://github.com/openvax/isovar/blob/master/isovar/variant_sequence.py):
Overlapping `AlleleRead`s containing the same mutation are assembled into a longer
sequence. The `VariantSequence` object represents this candidate coding sequence, as well
as all the `AlleleRead` objects which were used to create it.

* [ReferenceContext](https://github.com/openvax/isovar/blob/master/isovar/reference_context.py): To determine the reading frame in which to translate a `VariantSequence`, Isovar
looks at all Ensembl annotated transcripts overlapping the locus and collapses them
 into one or more `ReferenceContext` object. Each `ReferenceContext` represents the 
 cDNA sequence upstream of the variant locus and in which of the {0, +1, +2} reading frames
  it is translated. 

* [Translation](https://github.com/openvax/isovar/blob/master/isovar/translation.py): Use the reading frame from a `ReferenceContext` to translate a `VariantSequence` 
into a protein fragment, represented by `Translation`.

* [ProteinSequence](https://github.com/openvax/isovar/blob/master/isovar/protein_sequence.py):
Multiple distinct variant sequences and reference contexts can generate the same translations, so we aggregate those equivalent `Translation` objects into a `ProteinSequence`.

* [IsovarResult](https://github.com/openvax/isovar/blob/master/isovar/isovar_result.py): Since a single variant locus might have reads which assemble into multiple incompatible coding sequences, an `IsovarResult` represents a variant and one or more `ProteinSequence` objects which are associated with it. We typically don't want to deal with *every* possible translation of *every* distinct sequence detected around a variant, so the protein sequences are sorted by their number of supporting fragments and the best protein sequence is made easy to access. The `IsovarResult` object also has many informative properties such `num_alt_fragments`, `fraction_ref_reads`, &c.  



## Sequencing Recommendations

Isovar works best with high quality / high coverage mRNA sequence data. 
This means that you will get best results from >100M paired-end reads sequenced on an 
Illumina HiSeq from a library enriched with poly-A capture. The number of reads varies 
depending on degree of RNA degradation and tumor purity. The read length will determine 
the longest protein sequence you can recover, since Isovar's cDNA assembly only 
considers reads that overlap a variant. With 100bp reads you will be able to assemble
at most 199bp of sequence around a somatic single nucleotide variant, and consequently 
only be to determine 66 amino acids from the protein sequence. If you disable the cDNA 
assembly algorithm then a 100bp read will only be able to determine 33 amino acids.

## Other Commandline Tools

* `isovar-protein-sequences --vcf variants.vcf --bam rna.bam`

All protein sequences which can be assembled from RNA reads for any of the given variants.

* `isovar-allele-counts --vcf variants.vcf --bam rna.bam`

Counts of reads and fragments supporting the *ref*, *alt*, and *other* alleles at all given variant locations.

* `isovar-allele-reads --vcf variants.vcf --bam rna.bam`

Sequences of all reads overlapping any of the given variants. 
 
* `isovar-translations --vcf variants.vcf --bam rna.bam`

All possible translations of any assembled cDNA sequence containing any of the given variants in the reference frame of any matching transcript.
* `isovar-reference-contexts --vcf variants.vcf`

Shows all candidate reference contexts (sequence and reading frame) before each variant, derived from overlapping reference coding transcripts. 

* `isovar-variant-reads --vcf variants.vcf --bam rna.bam`

Like `isovar-allele-reads` but limited only to reads which support the *alt* allele.

* `isovar-variant-sequences --vcf variants.vcf --bam rna.bam`

Shows all assembled cDNA coding sequences supporting any of the given variants. 
