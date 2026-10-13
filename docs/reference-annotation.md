# Choosing the reference annotation

Isovar takes transcript models, reading frames and gene names from a
[PyEnsembl](https://github.com/openvax/pyensembl) annotation. Install the
Ensembl release that matches your alignments before the first run, for example:

```sh
pyensembl install --release 75 --species human
```

Align the RNA to the same assembly as the annotation.

## On the command line

`--genome` names an assembly such as `GRCh38` or `hg19`, and the most recent
Ensembl release installed for it is used. To choose the release, add it after
the assembly: `--genome GRCh38:93` or `GRCh38.93` (as in Ensembl's file names)
uses Ensembl 93. A release that does not provide the assembly, such as
`GRCh38:75`, is an error.

## From Python

Pass the same name (`genome="GRCh38:93"`) or a PyEnsembl genome when loading
the variants:

```python
import pyensembl
import varcode
from isovar import run_isovar

variants = varcode.load_vcf(
    "cancer-mutations.vcf", genome=pyensembl.EnsemblRelease(93))
isovar_results = run_isovar(variants=variants, alignment_file="tumor-rna.bam")
```

A `pyensembl.Genome` built from your own GTF and transcript FASTA files works
the same way, and gene and transcript names then come from that annotation.
