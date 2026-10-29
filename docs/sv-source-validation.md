# Validating SV source definitions

Inspect source rows before interpreting their RNA support:

```bash
isovar sv-source --vcf original.vcf --sha256 SOURCE_SHA256 \
  --reference-name GRCh38 --profile standard --output source-report.json.gz
```

`read_sv_source` is the corresponding typed Python API. Each `SVSourceRecord`
retains the original row and line number, explicit PCR state, caller identifiers,
unverified coordinate claims, interpretation reasons and any resolved Varcode
adjacency. `sv_source_report` creates the same JSON returned by the CLI. States
are `resolved_adjacency`, `unresolved` and `invalid`. They describe the input
definition, not RNA coverage, translation or tumor specificity.

`pcr_validation` is true, false or null. Null means no recorded PCR result. The
source filename and number of callers never change that value. Supplied sample
genotypes remain intact, including a benchmark nomination with genotype `0/0`.
An absent insertion sequence or endpoint is not assigned zero RNA support.

## Explicit source profiles

The default `standard` profile never repairs a malformed definition. Named
profiles require their documented original SHA-256 and GRCh38 assembly. They
must be selected explicitly; a similar filename or row is insufficient.

| Profile | Interpretation and remaining limitations |
|---|---|
| `colo829_paulin_paz_v2` | The four symbolic BND rows store literal oriented breakend alleles in END. Interpret those literal alleles through Varcode, require CHR2 agreement and record the transformation. Retain the zero-width inversion and seven sequence-unresolved insertions as unsupported definitions. |
| `hcc1395_svscope_2026` | Preserve 32 PCR-positive records separately from ensemble support. Keep missing/zero-width interval endpoints unresolved or invalid. Preserve caller-ID coordinates as unverified claims. The 22 added literal consensus-sequence records require reference anchoring; they are not ordinary positioned VCF alleles. |

The COLO829 source contains 54 rows: 46 resolved adjacency records, seven
unresolved records and one invalid inversion after the explicit interpretation.
The older 2022 inventory remains a separate source, including its own unresolved
insertions. Source rows are not necessarily independent biological events.

The HCC1395 source contains 1,118 rows: 157 resolved adjacency records, 780 invalid
interval definitions and 181 unresolved rows. None of its 22 added literal REF
sequences matches the local GRCh38 reference at the declared VCF position. Five
match 49 bases upstream; 17 have no exact match in the checked ±500-base region.
The published writer assigns a normal consensus to REF and a discovery-window
start to POS. Reference anchoring and any background differences must be resolved
before these sequences can be treated as positioned alleles.

Do not reconstruct missing endpoints as `POS + SVLEN`. For example, Minda_192
has `END=POS` and `SVLEN=11174710`, while two caller identifiers encode an endpoint
11,174,717 bases away. These claims remain separate, with no chosen repair.

## Annotation and RNA workflow

Source objects preserve their original contig spellings. Before predicting
annotation effects, use one resolver for the chosen annotation context:

```python
from pyensembl import EnsemblRelease
from varcode import ContigResolver
from isovar import read_sv_source

genome = EnsemblRelease(87)
resolver = ContigResolver(genome.contigs(), reference_name=genome.reference_name)
records = read_sv_source("original.vcf", sha256=source_digest,
                         reference_name="GRCh38", profile="standard")
for record in records:
    if record.variant is not None:
        annotation_variant = record.for_annotation(genome, resolver)
        effects = annotation_variant.effects()
```

Missing or ambiguous annotation contigs raise an explicit error. The original
source coordinates stay intact; aliases do not imply a liftover. Exact custom,
viral and decoy names are preserved during source parsing.

The [COLO829 manifest preparer](../reports/hcc1395_colo829/build_unified_sv_manifest.py)
uses this library interpretation and embeds it in each event's provenance. Its
resolved events enter the ordinary [SV RNA audit](sv-rna-audit.md); unresolved
rows retain explicit unassessed reasons. Source-specific interpretation does not
add a separate reconstruction or protein-inference path.

## Primary source conventions

- [VCF 4.5 specification](https://samtools.github.io/hts-specs/VCFv4.5.pdf):
  literal breakend orientation and symbolic interval fields.
- [COLO829 v2 source](https://zenodo.org/records/13917379), SHA-256
  `8ad075bf43aa71d89c591396569e813959ada44a30f171982c2d299411290bac`.
- [Minda amended COLO829 file](https://github.com/KolmogorovLab/minda/blob/6e13ee5120280d7a7fb6deef8c773a2d4cc41c8b/data/colo829_benchmark_grch38_amended_vaf.vcf):
  all four interpreted breakend definitions agree by coordinates and allele;
  event IDs differ across versions and are not used as correspondence evidence.
- [HCC1395 supplementary file](https://zenodo.org/records/18688406), SHA-256
  `892222606e3951f20a8608b28ff32e6dc320257d573c9b832108156e9a361424`.
- [Author's pre-PCR ensemble](https://github.com/Goatofmountain/SVScope/blob/0ad36860326490e94b8567276aa54535a79e5b26/analysis/minda_stat/Filter_2Tech5Tool/HCC1395_minda_ensemble.Filter.vcf):
  1,096 shared rows; interval END fields are absent in this version.
- [Author's VCF writer](https://github.com/Goatofmountain/SVScope/blob/0ad36860326490e94b8567276aa54535a79e5b26/src/OutVCF.py)
  and [SVScope paper](https://doi.org/10.1186/s13059-026-04076-0).

The [source-review receipts](../reports/sv_sources/README.md) include reproducible
reference comparisons and a replay of the four newly interpreted COLO829 events.

## Shared manifest preparation

`prepare_sv_source_record(record, genome, resolver, annotation="Ensembl 87")`
returns a typed `PreparedSVSourceRecord`: reference transcript models, every
Varcode hypothesis, both orientations of each adjacency, gene names, nominal
IGV loci and explicit annotation limitations. It performs no RNA reconstruction.
Its `SVSourceView.to_dict(reference_name)` and `reference_payload()` methods emit
the existing audit schemas. Reuse the same annotation resolver for the inventory.

Models and search views retain the original source contigs, including exact
custom names. A missing annotation contig leaves that RNA nomination available,
with no Varcode prediction and an explicit limitation. In contrast, the direct
`record.for_annotation` method raises when the requested annotation adaptation
cannot be completed. Neither path silently identifies an unknown contig as
intergenic. The helper selects all transcripts of genes in a declared breakpoint
window; reference completeness and observed assembly coverage remain separate.

[The general manifest preparer](../reports/sv_sources/build_sv_manifest.py) uses
this same helper for HCC1395 or another single source inventory. Supply a JSON
file with the audit schema's `files` and `sources` for the separate RNA libraries:

```bash
python reports/sv_sources/build_sv_manifest.py \
  --vcf original.vcf --sha256 SOURCE_SHA256 --reference-name GRCh38 \
  --profile hcc1395_svscope_2026 --sample-id HCC1395 \
  --libraries libraries.json --ensembl-release 87 \
  --region-flank 10000 --output prepared-inventory
isovar sv-audit --manifest prepared-inventory/manifest.json --output-dir rna-audit
```

The library file paths are relative to `libraries.json`; output manifests pin
absolute resolved paths and file digests. The report script only selects source
files and budgets and serializes library results. It retains every source row,
including unsupported definitions, and never pools cultures or processed libraries.
