# Isovar documentation

Isovar finds the mutant protein sequence that tumor RNA actually encodes around
each somatic mutation. This index lists the guides in reading order, from a
first run to specialized workflows. Start with
[Install](../README.md#install) and the [quick start](../README.md#quick-start)
in the project README.

## Getting started

| Guide | What it covers |
|---|---|
| [Reading and filtering results](results.md) | What each `IsovarResult` and table column means, and how filters mark results |
| [Command line](command-line.md) | Every `isovar` command, its output table, logging and errors |
| [Choosing the reference annotation](reference-annotation.md) | `--genome`, Ensembl releases and your own GTF |
| [Sequencing recommendations](sequencing-recommendations.md) | Library type, depth and read length, and what they limit |

## Running the pipeline

| Guide | What it covers |
|---|---|
| [Collecting and counting RNA reads](read-collection.md) | `ReadCollector` options, reads versus fragments, and the RNA support record in every JSON export |
| [Read eligibility across platforms](read-processing-audit.md) | Which alignments are used, and what Illumina, ONT and PacBio quality tags mean |
| [Shared RNA reconstruction and support](rna-reconstruction.md) | Public engine contracts across small variants, SVs, CTA isoforms and viral paths |
| [Assembly and translation](assembly-and-translation.md) | `ProteinSequenceCreator` options, and other variants found in the assembled RNA |
| [Protein context selection](protein-selection.md) | How much protein context is reported and which candidate comes first |
| [Reading frames from aligned reads](aligned-reading-frame.md) | How an indel upstream of the variant carries into the reading frame |
| [Phasing](phasing.md) | Which variants share RNA molecules, the cis/trans test, and matched germline variants |
| [Matched-normal evidence](matched-normal-evidence.md) | Reviewed normal inventories, genotype/reference-confidence states, supplied backgrounds and cohort integration |

## Exports and evidence for small variants

| Guide | What it covers |
|---|---|
| [Protein hypothesis export and Proteome adapter](protein-hypotheses.md) | Every protein and translation with scoped RNA evidence, including offline ORFdata composition |
| [RNA Proteomes](rna-proteome.md) | Shared API and CLI for SV, fusion and ncORF hypotheses, preserving unresolved alternatives |
| [Patient-background peptide pairs](patient-background-pairs.md) | Declared cDNA somatic-reversion baselines retaining germline changes, with separate reference-RNA support |
| [Competing allele interpretations](allele-interpretations.md) | Which of several callers' versions of one locus the RNA supports |
| [Mutation-evidence figures](visualization.md) | `isovar plot`, how to read its figures, and the osteosarc gallery |
| [Adapter and poly-A trimming](read-end-inference.md) | Opt-in annotation and trimming of technical read ends |
| [Base-quality evidence](base-quality.md) | Optional focal-allele quality categories reported alongside raw counts |
| [Partial-read support](partial-read-support.md) | Fractional support scores that divide each fragment among compatible hypotheses |
| [Cell/UMI evidence](cell-umi-evidence.md) | Cell barcode and UMI labels behind small-variant alleles and SV support |

## Structural variants and fusions

`run_isovar` accepts literal nucleotide alleles, including sequence-resolved
indels. Symbolic structural variants (`<DEL>`, `<DUP>`, etc.), breakends and
`varcode.StructuralVariant` objects are rejected rather than interpreted as
small variants. On the command line, SV records in a mixed `--vcf` are skipped
with a warning naming them. These workflows handle them:

| Guide | What it covers |
|---|---|
| [Validating SV source definitions](sv-source-validation.md) | `isovar sv-source`: original SV definitions and record-specific PCR evidence |
| [Unified SV RNA audit](sv-rna-audit.md) | `isovar sv-audit`: a source-separated manifest across nominated SVs, every RNA path assessed, competing hypotheses reconciled |
| [Reusing original RNA libraries](rna-library-index.md) | The library index dense event views reuse; `sv-audit --library-cache-dir` keeps verified indexes across runs |
| [SV RNA reconstruction](sv-rna.md) | `isovar sv-rna` / `reconstruct_sv_rna`: RNA paths around an SV call, their frames and ORFs, export and comparison |
| [Supplied fusion RNA](fusion.md) | `isovar fusion` / `reconstruct_fusion`: checking a fusion transcript assembled by another tool |
| [Literal-allele RNA paths](variant-rna-paths.md) | Intronic, coding and noncoding ALT/REF paths through the shared original-read and frame workflow |
| [ORF start evidence](orf-start-evidence.md) | Where an SV ORF's start codon comes from, in four tiers |
| [ONT read lineage](ont-read-lineage.md) | Split and duplex nanopore reads in SV support counts |
| [uORF evidence](uorf-evidence.md) | Offline human upstream-ORF models, protein sequences, native peptide evidence and mutation-position annotation |
| [Portable audit bundles](audit-bundles.md) | Installed Python API and CLI for preserving inputs, backup/restore and offline SV RNA replay |

## Viral sequences, antigen isoforms and reports

| Guide | What it covers |
|---|---|
| [Viral sequence recovery](viral-sequences.md) | `isovar viral`: observed viral RNA from GDC GRCh38.d1.vd1 alignments, without a mutation or human GTF |
| [CTA and other antigen isoforms](antigen-isoforms.md) | `isovar isoforms`: competing targeted transcripts, ambiguous RNA support and optional Salmon abundance |
| [Evidence report](../reports/hcc1395_colo829/PDF_REPORTS.md) | `isovar evidence-report`: a visual variant/RNA PDF and per-call TSV |
| [Antigen-source evidence report](antigen-evidence-report.md) | `isovar antigen-report`: SV, fusion, viral and antigen-isoform exports combined, with separate evidence units |
| [Reproducible RNA cohort reports](cohort-reports.md) | `isovar cohort`: multi-library reanalysis and per-sample PDFs |

## Concepts and design

| Guide | What it covers |
|---|---|
| [Library responsibilities](library-responsibilities.md) | How Varcode, Isovar and Vaxrank divide the work |
| [Internal design](internal-design.md) | How it works, and the objects that carry reads to protein sequences |
| [RNA assessment before protein selection](variant-read-assessment.md) | Retain intronic and noncoding ALT RNA with explicit assembly and frame states |
| [RNA in supplied ORF frames](variant-rna-orfs.md) | Reconstruct ncORF and ordinary coding candidates with original RNA support |
| [Observed-path reading frames](transcript-path-frames.md) | Shared transcript anchors, retained-intron interpretations and allele-versus-stop outcomes |

## Development and project records

- [Changelog](../CHANGELOG.md): behavior changes by release.
- [Releasing](../RELEASING.md): how a version is published.
- [Historical experiments](archive/experiments/README.md): inputs, measurements
  and code decisions behind earlier changes.
- The [notes](../notes) directory holds maintainer working notes: Sid analysis
  reports, design specifications and benchmark records, with their data. They
  are not usage guides.
- [Sid test reads](../notes/sid-test-reads.md): the original Sid osteosarcoma
  reads packaged with Isovar for its tests, and how to make them again.
- [Tumor-RNA audit](../tests/data/osteosarc/SAMPLE_AUDIT.md): a
  sample-by-sample RNA and protein audit of six osteosarc loci in two RNA
  sources, made for Isovar 1.8.0.
- [WLOS_0001 canine TP53 example](../tests/data/wlos0001_tp53/README.md): a
  published matched-normal tumor DNA call replayed against original,
  region-extracted tumor RNA reads in a small offline fixture.
- [HCC1395 and COLO829 source audit](../reports/hcc1395_colo829/README.md):
  published somatic call inventories, including SVs, and COLO829 polished
  isoform counts, with explicit RNA-support status for every source row.

### Design specifications

Working specifications for in-progress or completed engineering work. They
record goals and constraints for maintainers and are not usage guides.

- [SV RNA and protein reconstruction plan](sv-protein-assembly-plan.md)
- [Unified nominated-SV RNA audit (#305, #541, #563)](unified-sv-audit-spec.md)
- [Exon anchored SV path completion (#555)](sv-exon-bridge-555-spec.md)
- [Dense exon-anchored SV discovery (#559)](sv-dense-exon-anchor-spec.md)
- [Stream complete SV templates in dense support (#562)](sv3-dense-template-stream-562-spec.md)
- [Mutation-supporting RNA beyond canonical coding contexts](noncanonical-rna-discovery-spec.md)
- [DDX17 RNA-path investigation and deletion-aware report alignment](ddx17-splice-and-deletion-report-spec.md)
- [Complete the two unresolved ONT views (#470)](ont-igkc-recovery-spec.md)
- [Recovery event report memory](ont-report-memory-spec.md)
- [Share exact interval searches within each support batch](ont-shared-query-spec.md)
- [Exact early rejection experiment](ont-support-refinement-spec.md)
- [Bound complete-result validation (#500)](ont-validation-retention-spec.md)
- [Stream audit JSON (#488)](stream-audit-json-spec.md)
- [Incremental audit result decoding (#496)](stream-audit-decoding-spec.md)
- [Stream CLI checkpoint verification (#498)](stream-cli-resume-spec.md)
- [Complete report ledgers as release assets (#501)](sid-report-release-assets-spec.md)
