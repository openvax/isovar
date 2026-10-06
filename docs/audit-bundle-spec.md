# Portable audit bundles (#440)

Provide an installed, documented Python API and `isovar audit-bundle` CLI.
Neither requires a repository checkout, an LLM session, or machine-specific
working paths. Preserve original scientific evidence and receipt bytes.

## Public contract

- `AuditBundle.create(destination, files, metadata=None)` freezes explicitly
  selected files under logical relative names. `open`, `path`, `verify`, `pack`
  and `unpack` provide ordinary Python operations and matching CLI commands.
- A versioned JSON manifest records SHA-256, byte counts, software identity,
  metadata and a content identity. Reject missing/changed files, invalid paths,
  incompatible schemas and existing output destinations. Never fetch inputs.
- `bundle_sid_audit(directory, destination)` imports the existing Sid inventory,
  frozen transcript references, receipts, original regional BAMs/indexes and
  checkpoint/reconstruction files. Follow explicit file pins, validate them,
  and preserve original absolute paths only as provenance. Retain all outcomes,
  including incomplete acquisition and unassessable events.
- Imported SV reconstruction plans use logical BAM/index/input/result names.
  `replay_sv_rna(bundle, run_id, output_directory)` calls the existing public
  reconstruction engine with frozen inputs/settings. Verify the entire bundle
  before computing, record current software, compare with the frozen result,
  and resume only a matching, intact checkpoint. A changed engine needs an
  explicit option and still must pass the result comparison.
- Replay preserves acquisition limitations and product/geometry/orientation
  identities. It does not reacquire missing inputs or infer scientific outcomes.
  Frozen transcript models suffice for SV reconstruction; additional annotation
  assets needed by other analyses can be bundled explicitly. This PR does not
  add a second biological consequence engine or a new small-variant pipeline.

## Validation and delivery

Test round-trip and archive relocation, manifest/input corruption, invalid
paths and incomplete archives, interrupted publication, software/settings
drift, resume corruption, incomplete/unassessable outcome preservation, and
real original-read SV replay after removing its source directory. Exercise the
installed wheel's Python API and CLI from outside the checkout with networking
disabled. Export and verify the existing focused Sid subsets in durable storage
and replay a real selected event against its frozen reconstruction. Run lint,
the full suite and GitHub CI; bump the patch version, merge and deploy.
