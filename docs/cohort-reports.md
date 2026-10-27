# Reproducible RNA cohort reports

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
`prior_assessments`) write `receipt.json` alongside their result. SV receipts
also pin the referenced per-path artifacts, and prior assessments are copied
into a recorded snapshot. Outputs from 1.62.0 without these additional-task
receipts cannot be trusted by filename alone: compute them into a new output
directory with the current version. Existing small-variant receipts remain
renderable if their declarations and result digests match. Missing analyses
remain `not_run`; an existing result with a missing or mismatched receipt is an
error, not an empty analysis.

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

Supply `germline_variants` as another pinned allele array to annotate observed
background differences. A nomination's `origin_conflict` records a supplied
normal-call overlap; such calls remain visible but are excluded from report
counts of tumor-specific changed contexts. This is not the complete
patient-background peptide-pair API proposed in #100.

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
- `isoforms`: `targets` names a pinned JSON array of `IsoformTarget` fields,
  including a nested `FusionReference`; supply `reference_identity` and any
  reconstruction `parameters`.
- `prior_assessments`: `input` pins an earlier assessment and `summary` explains
  its scope. Reports label it as earlier evidence, not a newly executed screen.

A source may have `unavailable_reason`. `pending_inputs` rows use `sample_id`,
`source`, and `state`; the report identifies these gaps rather than treating
missing data as negative evidence.
An additional task's `sample_id` is required and must match its source. Supplied
fusion provenance must also name that sample; a conflicting input is rejected.

## Evidence and outputs

Per-source outputs include SQLite checkpoints, every call in `variants.jsonl.gz`,
a per-call TSV, all retained protein candidates in FASTA, and a receipt with
input, parameter, software and output identities. Errors retain their traceback;
they never become zero-support observations. Resume refuses changed inputs,
parameters or scientific implementation. Completed exports use deterministic
gzip headers. Per-call checkpoints use WAL so a long-running batch can be read
without losing completed work.

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
Reference ncORF Ribo-seq/MS evidence does not demonstrate mutant translation in
the current sample. Local protein contexts are not automatically complete ORFs,
and RNA support does not establish antigen presentation.

At more than 750 ALT observations, the default task uses the existing direct
read-context route rather than quadratic overlap assembly. No reads are dropped;
the route is recorded per call. Set
`max_alt_observations_for_overlap_assembly` explicitly for a sensitivity pass and
keep that pass separate. This budget can limit the extent of inferred sequences.
