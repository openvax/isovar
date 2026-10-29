# Reproducible RNA cohort reports

Protein funnels distinguish all local contexts from verified focal amino-acid
changes. Reference-identical and unresolved comparisons remain visible. Reports
show contributor, exact coding-window, single-read and merged paired-template
counts separately, on an identified nucleotide path. See
[coding coverage and comparisons](protein-hypotheses.md#coding-coverage-and-sequence-comparisons)
for the shared API and TSV semantics.

`isovar cohort` runs the same reconstruction APIs across declared RNA sources,
then produces one PDF per sample. Each library, culture, timepoint and analysis
pass keeps its own denominator. A processing copy or alternative transcript is
never counted as independent replication. Figures and PDF support are included
in the normal installation.

```bash
isovar cohort --manifest cohort.json --output-dir results
isovar cohort --manifest cohort.json --output-dir results --resume
isovar cohort --manifest cohort.json --output-dir results --render-only
```

For larger runs, compute individual sources and render together:

```bash
isovar cohort --manifest cohort.json --output-dir results \
  --source-id tumor-direct-RNA --no-render
isovar cohort --manifest cohort.json --output-dir results \
  --source-id tumor-PacBio --no-render
isovar cohort --manifest cohort.json --output-dir results --render-only
```

These commands can run separately, but reconstruction may require substantial
memory. The manifest records the source paths and SHA-256 digests. Relative
paths resolve beside the original manifest. Preserve that manifest and its
inputs, together with the outputs. Rendering consumes verified output artifacts
and does not rerun reconstruction.

Both the Python `run_cohort(..., render_only=True)` API and `--render-only`
require the recorded sample, source/library, scientific task settings, reference
and input digests to match the manifest. Output bytes are checked against their
recorded receipts before a report is generated. Changed declarations fail before
the saved manifest or summary is replaced; they cannot relabel existing evidence.
Use a new output directory for a different scientific run. Input paths may move
without changing the digests, and task `title`, `summary`, `locus_annotations`
and source `title` are presentation fields. A presentation annotation file must
still match its declared digest when used.

Rendering small-variant results does not require the original BAM, BAM index or
variant inputs. Fusion rendering still needs the pinned supplied-fusion input for
its parental reference sequences; SV gene labels need the declared annotation
file. Rendering preserves the producer's recorded software identity even when
the installed reporting code is newer. `--resume` also checks the scientific
implementation before reusing completed analysis.

Additional task families (`fusions`, `sv_audits`, `viral`, `isoforms`, and
`variant_rna_paths`, and `prior_assessments`) write `receipt.json` alongside their result. SV receipts
also pin the referenced per-path artifacts, and prior assessments are copied
into a recorded snapshot. Outputs from 1.62.0 without these additional-task
receipts cannot be trusted by filename alone: compute them into a new output
directory with the current version. Existing small-variant receipts remain
renderable if their declarations and result digests match. Missing analyses
remain `not_run`; an existing result with a missing or mismatched receipt is an
error, not an empty analysis.

Small-variant rows with an absent, ambiguous or unavailable alignment contig
have `status: unresolved_chromosome` and structured `contig_resolution`
diagnostics. Their RNA counts remain unmeasured (blank in TSV), rather than
zero. Receipts count `unresolved_chromosomes` separately from unexpected
`errors`; either produces `completed_with_errors`. PDFs retain these inputs in
the nomination denominator and exclude them from RNA-depth distributions.
See [contig resolution](read-collection.md#contig-names-and-unmeasured-loci).
If a sample supplies a `source_ids` display subset, every listed source must be
declared and belong to that sample; the subset cannot reassign another sample's
evidence.

## Antigen reconstruction policy

`viral` and `isoforms` tasks use the same `read_collection` mapping as the other RNA paths, including checksum-pinned read-end profiles. Their defaults retain MAPQ-zero/secondary reads and missing qualities, while excluding duplicate-marked reads and soft clips. An explicit policy applies to both discovery and complete-input support. The runner reuses validated original-library indexes across compatible tasks; different source identities are never pooled. See [shared policy controls](antigen-isoforms.md#shared-read-policy-and-library-reuse). Isoform tasks also write `transcripts.fasta` and `proteins.fasta`, with both files pinned in their receipts.

## Minimal manifest

The format is `isovar.cohort_manifest.v1`. Replace each example digest with the
SHA-256 of the corresponding file. BAM indexes must be beside their BAMs and are
also pinned in the analysis receipt.

```json
{
  "schema": "isovar.cohort_manifest.v1",
  "description": "Original RNA; separate libraries and cultures",
  "samples": [{"sample_id": "example"}],
  "files": {
    "rna": {"path": "rna.bam", "sha256": "SHA256"},
    "calls": {"path": "calls.json", "sha256": "SHA256"},
    "orfs": {"path": "human.orfdata.json.gz", "sha256": "SHA256"}
  },
  "sources": [{
    "source_id": "tumor-direct-RNA",
    "sample_id": "example",
    "library_id": "original-run-accession",
    "bam": "rna",
    "scope": "Complete original genomic RNA alignment",
    "provenance": {"accession": "original-run-accession", "assembly": "GRCh38"}
  }],
  "small_variants": [{
    "source_id": "tumor-direct-RNA",
    "pass_id": "default",
    "variants": "calls",
    "nc_orf_catalog": "orfs",
    "ensembl_release": 87,
    "protein_settings": {"protein_sequence_length": 49},
    "read_collection": {"min_mapping_quality": 1}
  }]
}
```

`calls.json` is an array of literal genomic alleles, for example
`[{"contig":"14","start":92811641,"ref":"GAGA","alt":"G"}]`. Starts are
1-based, using the untrimmed supplied allele. Preserve VCF/filter/source
information in each row's `origins`; consolidate duplicate alleles before the
run. Nonliteral SV alleles belong in `sv_audits`, not this array.

Contig names pass through Varcode's existing UCSC-to-annotation conversion;
`chrM`, `M`, `chrMT` and `MT` resolve to the mitochondrial annotation while the
original spelling remains available for BAM lookup. Receipts record conversions
and the runner emits a warning. Unknown names remain explicit failures. An
explicit source assembly conflicting with the annotation is rejected: name
conversion is not coordinate liftover or proof of reference-sequence identity.
Absence from a normal variant-only callset does not establish a reference genotype.

Supply `normal_evidence` as the key of a pinned
[`isovar normal-evidence` receipt](matched-normal-evidence.md) to retain selected
and excluded normal sources, local reference-confidence states and accepted
background alleles. The older `germline_variants` allele array also remains
supported; choose one reconciled background input per task. A nomination's
`origin_conflict` records supplied normal ALT or conflicting evidence; these
calls remain visible but are excluded from tumor-specific changed-context
counts. Use the [protein hypothesis API](protein-hypotheses.md) for the separate
patient-background peptide comparisons.

For a fully pinned annotation, list pinned annotation files under `files` and
include their keys in the task's `reference_files`. The chosen PyEnsembl release
must already be installed. `nc_orf_catalog` accepts a species/assembly-compatible
Orfdata catalogue. Failed transcript bindings remain explicit; catalogue
coverage is not a claim that all published ncORFs have been tested.

## Other evidence families

Each additional task has `id`, `sample_id`, and `source_id`:

- `sv_audits`: `manifest` points to an existing
  `isovar.sv_rna_audit_manifest.v1` file. The shared event/exon assessment and
  hypothesis reconciliation run on every retained path and its witnesses.
  Optional `locus_annotations` maps event IDs to `genes` and `nominal_igv`.
- `fusions`: `input` points to a supplied-fusion JSON input containing the
  nucleotide path, parent references and original RNA observations.
  `parameters` are passed to `reconstruct_fusion`.
- `viral`: the source's BAM must already contain competitive viral placements;
  `genome_fasta` names the pinned reference, and `parameters` are passed to
  `reconstruct_viral_sequences`. A human-only or regional BAM does not perform
  whole-library viral discovery through this task.
- `isoforms`: select `gene_set: "oncoref-cta"` or an explicit `genes` list, and
  supply `reference` as in the example below. API, CLI and cohort runs use the
  same target preparation and preserve the same selection ledger. Alternatively,
  `targets` names a pinned JSON array of `IsoformTarget` fields, including a
  nested `FusionReference`; supply `reference_identity`. These two input modes
  cannot be mixed. Reconstruction settings belong in `parameters`.
- `prior_assessments`: `input` pins an earlier assessment and `summary` explains
  its scope. Reports label it as earlier evidence, not a newly executed screen.

```json
{
  "id": "sample-1-CTA",
  "sample_id": "sample-1",
  "source_id": "sample-1-RNA",
  "gene_set": "oncoref-cta",
  "reference": {
    "genome_fasta": "alignment-fasta",
    "gtf": "annotation-gtf",
    "annotation_name": "Ensembl",
    "annotation_version": 115,
    "assembly": "GRCh38",
    "alignment_reference": "custom"
  }
}
```

`genome_fasta` and `gtf` above are keys in the manifest's checksum-pinned file
registry. Annotation release and target selection enter receipt identity;
changing either requires a new result. The installed OncoRef version is also
part of the computation identity. See [CTA membership](antigen-isoforms.md#where-cta-membership-comes-from)
for the distinction between catalogue nomination and downstream admission.

A source may have `unavailable_reason`. `pending_inputs` rows use `sample_id`,
`source`, and `state`; the report identifies these gaps rather than treating
missing data as negative evidence.
An additional task's `sample_id` is required and must match its source. Supplied
fusion provenance must match the source's declared input sample; a conflicting
input is rejected.
If the report groups specimens under a patient or cell-line label, set the
source's `input_sample_id` explicitly to the supplied fusion's original sample
ID (for example `sample_id: "Sid"`, `input_sample_id: "T1"`). The mapping is part
of the recorded identity. It does not rename original RNA witnesses or combine
different libraries; changing it invalidates result reuse. Without this field,
the fusion sample must equal the source's `sample_id`.

## Evidence and outputs

Per-source outputs include SQLite checkpoints, every call in `variants.jsonl.gz`,
a per-call TSV, all retained protein candidates in FASTA, and a receipt with
input, parameter, software and output identities. Errors retain their traceback;
they never become zero-support observations. Resume refuses changed inputs,
parameters or scientific implementation. Completed exports use deterministic
gzip headers. Per-call checkpoints use WAL so a long-running batch can be read
without losing completed work.

Landscape support distributions count calls at or above each observed positive
ALT-fragment depth. Equal depths share one threshold, and points keep singleton
categories visible. Zero-support calls remain in category totals and do not enter
the logarithmic support distribution.

The PDF separates focal-allele support, context contributors, exact full-window
allele-observation spans, single-original-read spans and basewise coding depth.
An allele observation can merge overlapping mates; only `source_read_count == 1`
qualifies for the single-original-read metric. Counts concern the selected
nucleotide path and are not summed across alternative translations. Per-sample
`figures/<sample>/summary.json` includes the 20-span thresholds used in the PDF.

Orange marks altered focal residues and deleted-residue boundaries. Matching
reference residues remain neutral; supplied germline, co-somatic and
unclassified differences use distinct colors. Reference and same-RNA focal
counterfactual comparisons do not fill missing sequence from annotation.
Fusion parental comparisons follow the observed translation start and both
edges of any inserted junction sequence. Inserted residues have no invented
parental counterpart. A donor-anchored protein remains visible when no coding
acceptor model is available, with the acceptor comparison explicitly unavailable.
Split junction codons use only contributing parental bases for comparison. An
identified acceptor frame shift is labeled; changed downstream residues use the
focal-difference color, while identical residues remain neutral.
Reference identity includes both annotation and transcript ID. If an acceptor
frame's transcript ID cannot identify one coding reference across the supplied
annotations, the report marks that comparison ambiguous instead of selecting an
arbitrary sequence. Distinct donor annotation identities remain distinct in the
FASTA. Fusion artifacts live under `figures/<sample>/fusions/<task-id>/`, with
`coverage.png`, `coverage.svg`, `coverage.pdf` and `proteins.fasta`; repeated
events in different libraries cannot overwrite each other's figures or sequences.
The input table uses every recorded additional-task input, including RNA BAMs,
reference FASTAs, target models and declared reference files. Unavailable tasks
do not require nonexistent input descriptors to render their unavailable state.
Reference ncORF Ribo-seq/MS evidence does not demonstrate mutant translation in
the current sample. Local protein contexts are not automatically complete ORFs,
and RNA support does not establish antigen presentation.

At more than 750 ALT observations, the default task uses the existing direct
read-context route rather than quadratic overlap assembly. No reads are dropped;
the route is recorded per call. Set
`max_alt_observations_for_overlap_assembly` explicitly for a sensitivity pass and
keep that pass separate. This budget can limit the extent of inferred sequences.

For intronic or other literal alleles needing extended original-read context,
use [literal RNA path tasks](variant-rna-paths.md). They share the CLI/API input
and one original-library index per source and policy; conditional RNA frames
remain separate from verified mutant-protein funnels.
