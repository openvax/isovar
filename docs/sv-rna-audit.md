# Unified structural-variant RNA audit

`isovar sv-audit` runs the existing RNA reconstruction, assesses **every retained
path** against the original library, then reconciles local RNA evidence with
Varcode hypotheses, ordinary transcripts, and competing original alignments.
The same workflow accepts validated supplied-fusion inputs. Reports consume its
snapshots; they do not infer frames independently.

Dense views share an [original-template library index](rna-library-index.md)
automatically. Add `--library-cache-dir rna-indexes` to retain verified indexes
across runs. Cache location does not change source separation or evidence budgets.

```sh
python -m isovar sv-audit --manifest manifest.json --output-dir results
python -m isovar sv-audit --manifest manifest.json --output-dir results --resume
```

The public Python functions are `assess_sv_rna_paths(result, bam, references)`,
`hypotheses_from_varcode(effects)`, `reconcile_sv_rna_hypotheses(assessment,
hypotheses)`, and `run_sv_rna_audit(manifest_path, output_dir)`. Assessment accepts
original FASTQ ends through `raw_reads=(read_group, query_name, mate, sequence)`;
these add sequence evidence, never genomic placements. Thresholds and budgets
are recorded in every result.

## Evidence and interpretation

Each path has separate junction witnesses, local branch support, exon/CDS links,
exact collinear annotated anchors, annotated-frame proteins, exploratory ORFs,
and full-interval sequence support. Exon alignments are stored once; path links
refer to `exon_alignment_id` in `assessment.exon_alignments`. Same-end exon clips and opposite mates are
separate links. Reference models retain strand and assembly. Original fragment
identity is `(read_group, query_name)` within a source; copied alignments do not
add fragments. This is not a claim of independent biological molecules.

An exon-linked read matching shared downstream reference sequence cannot
transfer its frame to an upstream event branch. Only the existing library's
annotated-frame or supplied-fusion validation produces annotated protein
candidates. Completeness requires an annotated 5′ start and a terminal stop;
otherwise the candidate remains partial. Exploratory ATGs retain separate
initiation and translation uncertainty. A read-overlap chain proves local
sequence tiling, not a mature full-length isoform. Its `chain_search` reports
`found`, `no_chain_within_read_limit`, or `search_limit`, with the evaluated
state count and explicit budgets.

Reconciliation uses original reads that exactly cover each junction or local
branch window. Exact query matches are **sequence compatibility**, without a
new base-quality or mapping-quality confidence model. Core placement support,
MAPQ, flags, clip qualities and original provenance remain available separately.
Absence of an exact match can reflect errors, incomplete coverage or an
unexamined placement and does not prove absence.

| Evidence state | Meaning |
|---|---|
| `supported` | The observed local window matches and distinguishes the hypothesis from at least one informative alternative. |
| `contradicted_by_informative_evidence` | The hypothesis has both unique collinear exact flanks, but its intervening sequence differs from an observed window. |
| `mixed_support_and_contradiction` | Different observed windows support and conflict with the hypothesis; coexisting isoforms remain possible. |
| `indistinguishable_in_observed_RNA` | Observed sequence matches but does not discriminate among examined hypotheses. |
| `insufficient_RNA` | No informative original window establishes local compatibility or conflict. |
| `unassessed` | No materialized candidate sequence is available. |

All candidates survive reconciliation, including Varcode candidates without
cDNA. Predicted complete transcript loss keeps its empty cDNA and
`predicted_empty_transcript` state; local RNA cannot validate that absence. Local alignment alternatives carry original sequences, placements, MAPQ
and read identity; they do not imply a whole transcript or protein. A shared
window can match multiple transcripts. Reconciliation does not infer initiation,
observed translation, normal/tumor specificity, or somatic causation.

## Manifest and reproducibility

The schema is `isovar.sv_rna_audit_manifest.v1`. The committed
[COLO829 manifest](../reports/hcc1395_colo829/unified_sv_audit/manifest.json)
is a full example. Its sections are:

- `files`: identifiers mapped to relative paths, SHA-256, and source descriptions.
- `sources`: source, sample, library, culture, role and provenance; either a BAM
  or an explicit `unavailable_reason`. Optional original FASTQ declares its
  name convention and read group. Different processed copies of one library
  remain identifiable and must never be summed as independent evidence.
- `events`: source-defined IDs, complete provenance and source IDs; oriented
  `views` or an explicit `unassessed_reason` for unresolved definitions.
- Each view names reference-model and optional hypothesis files, explicit
  donor/acceptor/regions inputs, and optional reconstruction budget overrides.
  `supplied_fusion` names an existing `isovar fusion` input; sample, event,
  assembly and reference models must agree. Its validated evidence remains in
  the reconstruction result; additional BAM sequence support remains separate.
- `parameters`, `read_collection`, `assessment`: explicit common budgets and
  policies, using the existing library API arguments.

Every input file is verified before computation. The runner records installed
sequence-library versions, all isovar Python source digests, BAM-index digests,
and the exact manifest. Each immutable compressed event/source/view snapshot
has a receipt. Resume verifies the implementation, inputs, parameters and output
hash; changed runs require a new output directory. `--event-id` selects explicit
events and records that restriction in the report. `report.json` and
`report.tsv` include unassessed and no-candidate runs, not just positives.

The report fixtures are bounded extracts of original data. Their limits are
part of the evidence. They do not establish that unretained exons, alternative
placements, other libraries, or other cultures lack a candidate.

## Scientific references

The coordinate and alignment model follows the
[SAM specification](https://samtools.github.io/hts-specs/SAMv1.pdf) and
[VCF specification](https://samtools.github.io/hts-specs/VCFv4.5.pdf).
[Varcode](https://github.com/openvax/varcode) owns the VCF adjacency and
consequence semantics; the cohort preparation preserves original rows and
keeps unresolved definitions explicit. Competitive ordinary-transcript and
alternative-placement assessment follows the rationale in
[FusionInspector](https://pmc.ncbi.nlm.nih.gov/articles/PMC10261907/).
