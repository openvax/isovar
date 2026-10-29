# Matched-normal evidence

`isovar normal-evidence` assesses literal nominated alleles against an explicit
inventory of normal VCF/BCF products. It writes a new receipt containing local
genotype and reference-confidence evidence, supplied background alleles, and the
scope of every selected or excluded source. It does not infer tumor specificity
from a missing ALT record.

```bash
isovar normal-evidence --manifest normal.json --output normal-evidence.json.gz
```

The Python entry points are `isovar.assess_matched_normal`,
`isovar.NormalSource`, `isovar.NormalEvidencePolicy`, and
`isovar.run_normal_evidence_manifest`. The CLI and cohort workflow consume the
same assessment. `isovar.normal_sources.MatchedNormalEvidence` describes the
returned receipt; per-source and per-allele results have corresponding typed
schemas.

`isovar.normal_background_regions` merges focal neighborhoods with explicitly
supplied background territory for acquisition and assessment. Both paths resolve
this territory against each source header with the same Varcode resolver.

## Source identity and input manifest

The `isovar.normal_manifest.v1` manifest contains:

- `sample_id` and `reference_name`: the patient/cohort and explicit common
  assembly. This does not perform liftover or establish sequence identity.
- `variants`: a `path` and `sha256` for a JSON/JSON.gz array of literal alleles
  with `contig`, 1-based `start`, `ref`, and `alt`. Original nomination metadata
  is retained in the receipt. Duplicate normalized alleles are rejected.
- `sources`: `NormalSource` descriptors for the complete reviewed inventory,
  including excluded and unavailable products.
- `background_flank`: nonnegative bases queried on each side of each focal
  interval, default 0. This bounds the supplied background; it is not a complete
  normal proteome or whole-gene germline assessment.
- `additional_background_regions`: optional `path` and `sha256` for a JSON/JSON.gz
  array of `[contig, start, end]` intervals in zero-based, half-open coordinates.
  Supply distant exons or ncORF coding blocks when reconstructed transcripts need
  background beyond a local genomic neighborhood. These regions and focal flanks
  form the requested background territory; the receipt retains the merged regions,
  source resolution failures, clipping and acquired scope.
- `snapshot`: an optional public inventory/release identifier.

Paths resolve relative to the manifest. A source descriptor supplies:

| Field | Meaning |
|---|---|
| `source_id` | Stable public file/product identifier; unique in this inventory |
| `selected`, `selection_reason` | Explicit selection or exclusion, independent of filenames |
| `specimen_id`, `library_id` | Reviewed physical specimen and library assignment |
| `vcf_sample` | Exact normal VCF column name; no suffix or capitalization rule |
| `caller`, `processing_group`, `derived_from` | Calling product and shared input/processing relationships |
| `reference_name` | Reviewed assembly declaration |
| `evidence_kind` | `variant_only` or `reference_confidence` |
| `path`, `sha256`, `index_path`, `index_sha256` | Local indexed VCF/BCF and both content digests |
| `acquired_regions` | Actual source-contig names with 0-based, half-open intervals; `null` means the complete file, `[]` means no territory |
| `policy` | Explicit genotype quality, depth and filtering requirements |
| `provenance` | Original URL, source headers, acquisition identity and assignment evidence |

Selected sources require all specimen/library/sample/caller assignments and
pinned file/index identities. Excluded sources can retain unresolved assignments.
Regional derivatives must declare their acquired intervals even if they retain
long records extending outside those intervals. Geometry is indexed once per
source; Varcode's contig resolver is also reused per header. Exact names,
mitochondrial aliases, alt and decoy contigs follow the same resolver rules as
[RNA read collection](read-collection.md#contig-names-and-unmeasured-loci).

For example, after constructing reviewed descriptors:

```python
from isovar import NormalEvidencePolicy, NormalSource, assess_matched_normal
from varcode import Variant

policy = NormalEvidencePolicy(minimum_genotype_quality=20, minimum_depth=5)
# source_descriptors is a reviewed list of the fields above, including digests.
sources = [NormalSource(**{**descriptor, "policy": policy}) for descriptor in source_descriptors]
receipt = assess_matched_normal(
    [Variant("chr1", 105, "A", "G", genome="GRCh38")],
    sources, sample_id="patient", reference_name="GRCh38", background_flank=500,
    additional_background_regions=[("chr1", 800, 900)])
```

## Evidence states and quality policy

| State | Conclusion permitted by the supplied evidence |
|---|---|
| `alternate_observed` | An accepted normal genotype calls this exact literal allele |
| `conflicting_calls` | Accepted reference and alternate evidence disagree within or across products |
| `confidently_reference` | Accepted reference genotypes cover every tested base |
| `other_allele_observed` | An accepted non-reference genotype overlaps this interval, but does not call this exact allele |
| `uncalled` | Queried records are unqualified/no-call, or a reference-confidence source has no covering call |
| `variant_only_absence` | No overlapping record in an acquired variant-only source |
| `unavailable` | No selected evidence, missing/ambiguous source contig, unacquired territory or a failed query |

The reconciled state always accompanies its per-source assessments. An ALT call
in one product does not erase unavailable evidence in another. Excluded products
and their reasons are stored once in the source inventory, rather than repeated
for every allele. Processing copies never vote or count as independent specimens.

The default analysis policy requires PASS record/sample filters, a complete
genotype, GQ >= 20 and DP >= 5. For callers such as Strelka, explicitly choose
`genotype_quality_field: "GQX"` and `indel_depth_field: "DPI"`. Strelka's DPI is
depth at the position preceding an indel, not per-base interval coverage; it is
used only for an explicitly declared indel-genotype gate. Reference blocks
always require `MIN_DP`; average block
DP is insufficient. Unfiltered `FILTER=.` is rejected unless
`allow_unfiltered: true` is explicitly supplied. These thresholds are analysis
choices, not universal biological cutoffs. Missing confidence fields remain
unresolved; they are never silently replaced by site QUAL.

Reference confidence covers the changed reference bases, or both flanks of an
insertion. Adjacent reference blocks can cover an interval; gaps cannot. A
symbolic SV's 0/0 genotype does not establish base-level reference over its END
interval. Explicit overlapping REF disagreement invalidates that record.

Only genotype-indexed, accepted literal ALT alleles enter `background_variants`.
Background collection makes one sequential pass per selected source, independently
of focal query success. A complete observed record is selected when its literal
REF span overlaps the intersection of requested and acquired territory. This
preserves boundary indels without claiming reference confidence outside that
territory. A symbolic SV's END does not expand the literal REF span for selection.
Focal reference-confidence assessment still requires coverage of every tested base.
For example, `REF=A, ALT=G,T, GT=0/2` supplies T, not G. Original genotype,
phase set, allelic depths, filters, confidence, caller/source and processing
provenance remain available. Supplied normal phasing is not evidence that an RNA
assembly traversed those alleles. Alleles are trimmed with Varcode; repeat-shifted
indels are not automatically left-aligned or treated as equivalent haplotypes.

## Cohort integration and versioning

Add the receipt as a pinned cohort `files` entry and set the small-variant task's
`normal_evidence` to that entry's key. It supplies background alleles and attaches
the per-nomination normal assessment. A missing nomination in that receipt is
explicitly unavailable. Accepted focal ALT or conflicting calls are flagged for
origin review and excluded from tumor-specific changed-context funnels.

The older `germline_variants` array remains supported, but cannot be combined
with `normal_evidence`: reconcile those inputs first. Existing nomination origin
flags are preserved. To reassess an earlier normal union, start from the original
unflagged nominations and write a new receipt and a new cohort output directory.
Changed normal evidence invalidates resume and render identity; completed
scientific receipts are not rewritten.

Cohort PDFs show normal-state counts, selected specimens/libraries/callers and
policies, excluded-product counts/reasons, background scope and the receipt
digest. Full source IDs and record provenance remain in the pinned normal
receipt. The RNA cohort receipt embeds its source inventory for offline rendering.

The [Sid inventory](../reports/matched_normal/sid_normal_inventory.json.gz) pins
the reviewed public source selection. Its acquisition script downloads original
public file streams and filters locally, hashing the complete stream and indexed
derivative. It sends no patient coordinates in HTTP queries. See the adjacent
report README for the reproducible analysis and its actual limitations.

Primary definitions: [VCF specification](https://samtools.github.io/hts-specs/VCFv4.3.pdf)
and [GATK reference-confidence blocks](https://gatk.broadinstitute.org/hc/en-us/articles/360035531812-GVCF-Genomic-Variant-Call-Format).
Strelka's caller-specific fields are documented by
[Illumina](https://support.illumina.com/help/BaseSpace_App_WGS_v5_OLH_15050955_02/Content/Source/Informatics/Apps/VCFFileAnnotations_appENR_appISCWGS.htm).
