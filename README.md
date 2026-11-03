[![Tests](https://github.com/openvax/isovar/actions/workflows/tests.yml/badge.svg)](https://github.com/openvax/isovar/actions/workflows/tests.yml)
[![Coverage Status](https://coveralls.io/repos/openvax/isovar/badge.svg?branch=master&service=github)](https://coveralls.io/github/openvax/isovar?branch=master)
[![PyPI](https://img.shields.io/pypi/v/isovar.svg?maxAge=1000)](https://pypi.org/project/isovar/)

# Isovar

Isovar finds the mutant protein sequence that tumor RNA actually encodes around
each somatic mutation. Given variant calls (VCF) and aligned tumor RNA-seq reads
(BAM), it assembles the reads that carry each mutation into cDNA and translates
that cDNA in the reading frames of annotated transcripts.

The protein comes from the reads, not from the reference plus one edit, so it
can include nearby variants and splice junctions that the reads show. Where
coverage runs out or no reading frame can be established, the result is shorter
or empty: Isovar never fills gaps with reference sequence.

Isovar is one of three OpenVax libraries that share this work.
[Varcode](https://github.com/openvax/varcode) generates transcript hypotheses
and predicts coding consequences; Isovar reconstructs RNA-supported sequences
and reconciles the evidence; [Vaxrank](https://github.com/openvax/vaxrank)
evaluates protein and peptide candidates. See
[library responsibilities](https://github.com/openvax/isovar/blob/master/docs/library-responsibilities.md).

## Install

```sh
pip install isovar
```

<a id="installation"></a>Isovar requires Python 3.11 or later. Reference annotation comes from
[PyEnsembl](https://github.com/openvax/pyensembl); install the Ensembl release
that matches your alignments before the first run, for example:

```sh
pyensembl install --release 75 --species human
```

To choose a specific release or use your own GTF, see
[choosing the reference annotation](https://github.com/openvax/isovar/blob/master/docs/reference-annotation.md).
Figures and PDF reports are included in the standard installation.

## Quick start

```python
from isovar import run_isovar, isovar_results_to_dataframe

isovar_results = run_isovar(
    variants="cancer-mutations.vcf",
    alignment_file="tumor-rna.bam")

for isovar_result in isovar_results:
    # The protein preferred by the context/support policy, or None.
    protein = isovar_result.top_protein_sequence
    if protein is not None:
        # num_alt_fragments: distinct fragments supporting the variant allele.
        print(isovar_result.variant, isovar_result.num_alt_fragments, protein.amino_acids)

# One row per variant.
df = isovar_results_to_dataframe(isovar_results)
```

With this repository's B16 mouse melanoma test files
(`tests/data/b16.f10/b16.vcf` and `b16.combined.sorted.bam`, after
`pyensembl install --release 102 --species mouse`), it prints:

```text
Variant(contig='9', start=82927102, ref='G', alt='T', reference_name='GRCm38') 15 GIEDKITKMKVTMVAWDRHDNTVIHAVNNMTLKVWNSYTGQLIHVLMGH
Variant(contig='X', start=8125624, ref='C', alt='A', reference_name='GRCm38') 21 FDVERAAKAAVNKLQGHSAPVLDVIVNCDESLLASSDASGMVIVWRREQ
```

The file's other two variants have no supporting RNA fragments and no protein.
They are still in `isovar_results` and in `df`, where they fail the default
filters.

<a id="reading-the-results"></a>
<a id="filtering-results"></a>

Each result holds the RNA evidence at one variant's locus and the proteins
assembled from it. Filters mark results but never remove them. See
[reading and filtering results](https://github.com/openvax/isovar/blob/master/docs/results.md)
for every property, table column and filter option.

Isovar logs through the standard `logging` module under the `isovar` logger and
never configures logging itself; configure it in your application to see
progress.

## Command line

```sh
isovar run \
    --vcf somatic-variants.vcf \
    --bam rnaseq.bam \
    --output isovar-results.csv
```

The CSV has one row per variant, with the same columns as the DataFrame.
`isovar --help` lists the other commands: table exports, RNA hypotheses,
structural variants and fusions, viral and antigen-isoform evidence, proteomes,
reports and figures, and audits. The
[command-line guide](https://github.com/openvax/isovar/blob/master/docs/command-line.md)
describes each one and its output.

## Documentation

<a id="project-reports-and-test-data"></a>The [documentation index](https://github.com/openvax/isovar/blob/master/docs/README.md)
lists every guide in reading order. Start with:

- [Reading and filtering results](https://github.com/openvax/isovar/blob/master/docs/results.md)
- <a id="sequencing-recommendations"></a>[Choosing the reference annotation](https://github.com/openvax/isovar/blob/master/docs/reference-annotation.md)
  and [sequencing recommendations](https://github.com/openvax/isovar/blob/master/docs/sequencing-recommendations.md)
- [Collecting and counting RNA reads](https://github.com/openvax/isovar/blob/master/docs/read-collection.md)
  and [protein context selection](https://github.com/openvax/isovar/blob/master/docs/protein-selection.md)
- <a id="structural-variants-and-fusions"></a>Structural variants: `run_isovar` rejects symbolic SVs and breakends; use
  [`isovar sv-rna`](https://github.com/openvax/isovar/blob/master/docs/sv-rna.md),
  [`isovar sv-audit`](https://github.com/openvax/isovar/blob/master/docs/sv-rna-audit.md)
  or, for a fusion assembled by another tool,
  [`isovar fusion`](https://github.com/openvax/isovar/blob/master/docs/fusion.md)
- <a id="how-it-works"></a>[Internal design](https://github.com/openvax/isovar/blob/master/docs/internal-design.md):
  how reads become proteins
- [Changelog](https://github.com/openvax/isovar/blob/master/CHANGELOG.md)

## Development

```sh
./develop.sh   # editable install with uv
./lint.sh      # ruff
./test.sh      # pytest
```

`./test.sh` needs the Ensembl releases that
[`.github/scripts/install_ensembl_data.sh`](https://github.com/openvax/isovar/blob/master/.github/scripts/install_ensembl_data.sh)
installs (human 75, 87 and 95, mouse 102). See
[RELEASING.md](https://github.com/openvax/isovar/blob/master/RELEASING.md) for
version bumps and releases.

## Citation

<a id="citing-isovar"></a>Isovar is part of the OpenVax pipeline described in:

> Rubinsteyn, A., Kodysh, J., Hodes, I., et al.
> [Computational Pipeline for the PGV-001 Neoantigen Vaccine Trial.](https://www.frontiersin.org/journals/immunology/articles/10.3389/fimmu.2017.01807/full)
> *Frontiers in Immunology* 8, 1807 (2018).

> Kodysh, J. & Rubinsteyn, A.
> [OpenVax: An Open-Source Computational Pipeline for Cancer Neoantigen Prediction.](https://link.springer.com/protocol/10.1007/978-1-0716-0327-7_10)
> *Methods Mol. Biol.* 2120, 147–160 (2020).

## License

Isovar is released under the [Apache License 2.0](https://github.com/openvax/isovar/blob/master/LICENSE).
