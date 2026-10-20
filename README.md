[![Tests](https://github.com/openvax/isovar/actions/workflows/tests.yml/badge.svg)](https://github.com/openvax/isovar/actions/workflows/tests.yml)
<a href="https://coveralls.io/github/openvax/isovar?branch=master">
    <img src="https://coveralls.io/repos/openvax/isovar/badge.svg?branch=master&service=github" alt="Coverage Status" />
</a>
<a href="https://pypi.python.org/pypi/isovar/">
    <img src="https://img.shields.io/pypi/v/isovar.svg?maxAge=1000" alt="PyPI" />
</a>

# Isovar

Isovar finds the mutant protein sequence that tumor RNA actually encodes around
each somatic mutation. Given variant calls (VCF) and aligned tumor RNA-seq reads
(BAM), it assembles the reads that carry each mutation into cDNA and translates
that cDNA in the reading frames of annotated transcripts.

The protein comes from the reads, not from the reference plus one edit, so it
can include nearby variants and splice junctions that the reads show. Where
coverage runs out or no reading frame can be established, the result is shorter
or empty: Isovar never fills gaps with reference sequence.

## Installation

```sh
pip install isovar
```

Isovar requires Python 3.11 or later. Reference annotation comes from
[PyEnsembl](https://github.com/openvax/pyensembl); install the Ensembl release
that matches your alignments before the first run, for example:

```sh
pyensembl install --release 75 --species human
```

To choose a specific release or use your own GTF, see
[choosing the reference annotation](https://github.com/openvax/isovar/blob/master/docs/reference-annotation.md).
The figure commands (`isovar plot` and `isovar fusion --plot-dir`) need the
plotting extra: `pip install 'isovar[plot]'`.

## Quick start

From Python:

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

Isovar logs through the standard `logging` module under the `isovar` logger and
never configures logging itself; configure it in your application to see
progress.

From the command line:

```sh
isovar run \
    --vcf somatic-variants.vcf \
    --bam rnaseq.bam \
    --output isovar-results.csv
```

The CSV has one row per variant, with the same columns as the DataFrame.
`isovar --help` lists the other commands; the
[command-line guide](https://github.com/openvax/isovar/blob/master/docs/command-line.md)
describes each one and its output.

## Reading the results

`run_isovar` returns one `IsovarResult` per input variant, in input order. Each
holds the RNA evidence at that variant's locus and any mutant protein sequences
assembled for it.

| `IsovarResult` property | Table column | Meaning |
|---|---|---|
| `variant` | `variant`, `chr`, `pos`, `ref`, `alt` | The input variant |
| `top_protein_sequence` | `protein_sequence` and other `protein_sequence_*` columns | The top-ranked mutant protein assembled from RNA, or none |
| `sorted_protein_sequences` | | Every protein kept, ranked; one by default |
| `num_alt_reads`, `num_alt_fragments`, `fraction_ref_reads`, ... | Same names | Read evidence for the ref, alt and other alleles |
| `predicted_effect` | `predicted_effect` | Varcode's prediction from the reference plus this one edit, for comparison |
| `passes_all_filters`, `filter_values` | `passes_all_filters`, `filter:<name>` | Filter outcomes |

By default Isovar aims for 49 residues of protein: 2K−1 for peptide size
K = 25, enough to contain every 25-mer overlapping a single changed residue.
The actual length depends on RNA coverage.
[Protein context selection](https://github.com/openvax/isovar/blob/master/docs/protein-selection.md)
explains the target and how the top protein is chosen.

A read is one sequenced segment, such as one mate of a pair; a fragment is one
template, so both mates of a pair count once. Fragment counts are usually the
ones to filter on. Neither is a molecule count; see
[collecting and counting RNA reads](https://github.com/openvax/isovar/blob/master/docs/read-collection.md).

### Filtering results

`run_isovar` evaluates filters on each result; a failing result is kept, with
`False` in its `filter_values` dictionary and in `passes_all_filters`. In the
DataFrame and CSV each filter becomes a `filter:<name>` column. Filters never
remove a result or a row.

`filter_thresholds` maps names like `'min_num_alt_reads'` or
`'max_fraction_other_fragments'` to numbers. The text after `min_` or `max_`
names a numeric property of `IsovarResult`, and most read-evidence properties
follow the pattern `{num|fraction}_{ref|alt|other}_{reads|fragments}`. For
example, this requires at least 10 alt reads and at most 25% of fragments
supporting other alleles:

```python
from isovar import run_isovar

isovar_results = run_isovar(
    variants="cancer-mutations.vcf",
    alignment_file="tumor-rna.bam",
    filter_thresholds={"min_num_alt_reads": 10, "max_fraction_other_fragments": 0.25})

for isovar_result in isovar_results:
    print(isovar_result.variant, isovar_result.passes_all_filters)
```

`filter_flags` names boolean properties of `IsovarResult`; prefix one with
`not_` to negate it, as in `not_protein_sequence_matches_predicted_mutation_effect`.

Omitting either argument applies the defaults in
[`default_parameters.py`](https://github.com/openvax/isovar/blob/master/isovar/default_parameters.py):
`DEFAULT_FILTER_THRESHOLDS`, and `DEFAULT_FILTER_FLAGS`, which are
`predicted_effect_modifies_protein_sequence`,
`has_mutant_protein_sequence_from_rna` and `protein_sequence_contains_mutation`.
Passing a value replaces the corresponding defaults; to change one threshold,
copy `DEFAULT_FILTER_THRESHOLDS` and update it.

## Structural variants and fusions

The pipeline above accepts literal nucleotide alleles, including
sequence-resolved indels. Symbolic structural variants (`<DEL>`, `<DUP>`, etc.),
breakends and `varcode.StructuralVariant` objects are rejected rather than
interpreted as small variants. Two separate workflows handle them:

- `isovar sv-rna` / `reconstruct_sv_rna` reconstructs exploratory RNA paths
  around one nominated SV from a BAM and annotated transcript models
  ([guide](https://github.com/openvax/isovar/blob/master/docs/sv-rna.md)).
- `isovar fusion` / `reconstruct_fusion` checks the junction evidence and
  annotated coding frames of a fusion transcript assembled by another tool
  ([guide](https://github.com/openvax/isovar/blob/master/docs/fusion.md)).

## How it works

For antigen genes without a nominated mutation, `isovar isoforms` reconstructs
competing CTA or other targeted transcripts, retains ambiguous RNA support,
and optionally estimates abundance with Salmon. See
[antigen isoforms](https://github.com/openvax/isovar/blob/master/docs/antigen-isoforms.md)
and the [preferred combined human/viral reference](https://github.com/openvax/isovar/blob/master/docs/reference-annotation.md#combined-human-and-viral-reference).

For RNA aligned to GDC GRCh38.d1.vd1, `isovar viral` recovers observed viral
sequences across its 200-contig panel without a mutation or human GTF. It
exports partial/spliced nucleotide paths and original fragment evidence.
Sequence-qualified ORFdata annotations also provide predicted coding products
and protein FASTA, initially covering nine RefSeq HPV16 CDSs. Sample variants
and partial coding ends remain explicit. See
[viral sequence recovery](https://github.com/openvax/isovar/blob/master/docs/viral-sequences.md).

For each variant, Isovar:

1. collects the RNA reads that overlap it,
2. keeps the reads that carry the mutant allele,
3. assembles overlapping mutant reads into longer cDNA sequences,
4. places those sequences in the reading frames of annotated transcripts, and
5. translates them into mutant protein sequences.

[Internal design](https://github.com/openvax/isovar/blob/master/docs/internal-design.md)
describes the objects behind each step.

Isovar is one of three OpenVax libraries that share this work.
[Varcode](https://github.com/openvax/varcode) generates transcript hypotheses
and predicts coding consequences; Isovar reconstructs RNA-supported sequences
and reconciles the evidence; [Vaxrank](https://github.com/openvax/vaxrank)
evaluates protein and peptide candidates. See
[library responsibilities](https://github.com/openvax/isovar/blob/master/docs/library-responsibilities.md).

## Documentation

**Running the pipeline**

| Guide | What it covers |
|---|---|
| [Choosing the reference annotation](https://github.com/openvax/isovar/blob/master/docs/reference-annotation.md) | `--genome`, Ensembl releases and your own GTF |
| [Command line](https://github.com/openvax/isovar/blob/master/docs/command-line.md) | Every `isovar` command, its output table, logging and errors |
| [Collecting and counting RNA reads](https://github.com/openvax/isovar/blob/master/docs/read-collection.md) | `ReadCollector` options, reads versus fragments, and the RNA support record in every JSON export |
| [Read eligibility across platforms](https://github.com/openvax/isovar/blob/master/docs/read-processing-audit.md) | Which alignments are used, and what Illumina, ONT and PacBio quality tags mean |
| [Assembly and translation](https://github.com/openvax/isovar/blob/master/docs/assembly-and-translation.md) | `ProteinSequenceCreator` options, and other variants found in the assembled RNA |
| [Protein context selection](https://github.com/openvax/isovar/blob/master/docs/protein-selection.md) | How much protein context is reported and which candidate comes first |
| [Reading frames from aligned reads](https://github.com/openvax/isovar/blob/master/docs/aligned-reading-frame.md) | How an indel upstream of the variant carries into the reading frame |
| [Phasing](https://github.com/openvax/isovar/blob/master/docs/phasing.md) | Which variants share RNA molecules, the cis/trans test, and matched germline variants |

**Exports and evidence for small variants**

| Guide | What it covers |
|---|---|
| [Protein hypothesis export](https://github.com/openvax/isovar/blob/master/docs/protein-hypotheses.md) | Every protein and translation with scoped RNA evidence, for tools that combine sources |
| [Competing allele interpretations](https://github.com/openvax/isovar/blob/master/docs/allele-interpretations.md) | Which of several callers' versions of one locus the RNA supports |
| [Mutation-evidence figures](https://github.com/openvax/isovar/blob/master/docs/visualization.md) | `isovar plot`, how to read its figures, and the osteosarc gallery |
| [Adapter and poly-A trimming](https://github.com/openvax/isovar/blob/master/docs/read-end-inference.md) | Opt-in annotation and trimming of technical read ends |
| [Base-quality evidence](https://github.com/openvax/isovar/blob/master/docs/base-quality.md) | Optional focal-allele quality categories reported alongside raw counts |
| [Partial-read support](https://github.com/openvax/isovar/blob/master/docs/partial-read-support.md) | Fractional support scores that divide each fragment among compatible hypotheses |
| [Cell/UMI evidence](https://github.com/openvax/isovar/blob/master/docs/cell-umi-evidence.md) | Cell barcode and UMI labels behind small-variant alleles and SV support |

**Structural variants and fusions**

| Guide | What it covers |
|---|---|
| [SV RNA reconstruction](https://github.com/openvax/isovar/blob/master/docs/sv-rna.md) | `isovar sv-rna`: RNA paths around an SV call, their frames and ORFs, export and comparison |
| [Supplied fusion RNA](https://github.com/openvax/isovar/blob/master/docs/fusion.md) | `isovar fusion`: checking a fusion transcript assembled by another tool |
| [ORF start evidence](https://github.com/openvax/isovar/blob/master/docs/orf-start-evidence.md) | Where an SV ORF's start codon comes from, in four tiers |
| [ONT read lineage](https://github.com/openvax/isovar/blob/master/docs/ont-read-lineage.md) | Split and duplex nanopore reads in SV support counts |
| [uORF evidence](https://github.com/openvax/isovar/blob/master/docs/uorf-evidence.md) | Offline human upstream-ORF models, protein sequences, native peptide evidence and mutation-position annotation |
| [Portable audit bundles](https://github.com/openvax/isovar/blob/master/docs/audit-bundles.md) | Installed Python API and CLI for preserving inputs, backup/restore and offline SV RNA replay |

**Project**

| Guide | What it covers |
|---|---|
| [Library responsibilities](https://github.com/openvax/isovar/blob/master/docs/library-responsibilities.md) | How Varcode, Isovar and Vaxrank divide the work |
| [Internal design](https://github.com/openvax/isovar/blob/master/docs/internal-design.md) | The objects that carry reads to protein sequences |
| [Changelog](https://github.com/openvax/isovar/blob/master/CHANGELOG.md) | Behavior changes by release |

## Sequencing recommendations

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

## Project reports and test data

The [notes](https://github.com/openvax/isovar/tree/master/notes) directory holds
maintainer working notes: Sid analysis reports, design specifications and
benchmark records, with their data. They are not usage guides.

- [Sid test reads](https://github.com/openvax/isovar/blob/master/notes/sid-test-reads.md):
  the original Sid osteosarcoma reads packaged with Isovar for its tests, and
  how to make them again.
- [Tumor-RNA audit](https://github.com/openvax/isovar/blob/master/tests/data/osteosarc/SAMPLE_AUDIT.md):
  a sample-by-sample RNA and protein audit of six osteosarc loci in two RNA
  sources, made for Isovar 1.8.0.

## Citing Isovar

Isovar is part of the OpenVax pipeline described in:

> Rubinsteyn, A., Kodysh, J., Hodes, I., et al.
> [Computational Pipeline for the PGV-001 Neoantigen Vaccine Trial.](https://www.frontiersin.org/journals/immunology/articles/10.3389/fimmu.2017.01807/full)
> *Frontiers in Immunology* 8, 1807 (2018).

> Kodysh, J. & Rubinsteyn, A.
> [OpenVax: An Open-Source Computational Pipeline for Cancer Neoantigen Prediction.](https://link.springer.com/protocol/10.1007/978-1-0716-0327-7_10)
> *Methods Mol. Biol.* 2120, 147–160 (2020).

## License

Isovar is released under the [Apache License 2.0](https://github.com/openvax/isovar/blob/master/LICENSE).
