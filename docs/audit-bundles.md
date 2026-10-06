# Portable audit inputs and offline replay

Install Isovar normally (`pip install isovar`). The public Python API and
`isovar audit-bundle` / `isovar-audit-bundle` commands work from any directory.
They do not need a checkout, a notebook, an LLM session or network access.

The [preserved Sid bundles and installed-wheel validation](sid-portable-audit-bundles.md)
provide a worked example with original sequencing inputs and recorded identities.

## Bundle explicitly selected files

```python
from isovar.audit_bundle import AuditBundle

bundle = AuditBundle.create(
    "saved-audit",
    files={
        "reads/tumor.bam": "regional-reads.bam",
        "reads/tumor.bam.bai": "regional-reads.bam.bai",
        "annotation/transcripts.json": "frozen-transcripts.json",
    },
    metadata={"sample": "tumor-T2", "purpose": "selected SV audit"},
)
print(bundle.verify())
print(bundle.path("reads/tumor.bam"))
bundle.pack("saved-audit.tar.gz")
restored = AuditBundle.unpack("saved-audit.tar.gz", "another-location")
assert restored.bundle_id == bundle.bundle_id
```

`create` copies only the named files. `expected_sha256={logical_name: checksum}`
can require existing input pins. Source files and destination directories must
be local. Existing destinations are never overwritten. Interrupted creation
does not publish a complete manifest. Keep the archive on persistent storage
and copy it to your backup destination using ordinary file-copy tools.

The equivalent CLI reads a JSON specification. File paths are relative to the
specification file, so it can be run from any working directory:

```json
{
  "files": {
    "reads/tumor.bam": "regional-reads.bam",
    "reads/tumor.bam.bai": "regional-reads.bam.bai"
  },
  "metadata": {"sample": "tumor-T2"}
}
```

```sh
isovar audit-bundle create --spec inputs.json --output saved-audit
isovar audit-bundle verify saved-audit
isovar audit-bundle pack saved-audit --output saved-audit.tar.gz
isovar audit-bundle unpack saved-audit.tar.gz --output another-location
```

## Import and replay an existing Sid audit

```sh
isovar audit-bundle import-sid --audit EXISTING_AUDIT --output saved-sid-audit
isovar audit-bundle inspect saved-sid-audit
isovar audit-bundle pack saved-sid-audit --output saved-sid-audit.tar.gz
isovar audit-bundle unpack saved-sid-audit.tar.gz --output relocated-sid-audit
isovar audit-bundle replay relocated-sid-audit --run-id SOURCE/EVENT/reverse \
  --output replay-checkpoints
```

`inspect` lists exact run IDs. Repeat `--run-id` for a subset or use `--all` for
every replayable view. Each source, event and orientation stays separate.
Completed unassessable/timeout outcomes are retained in `outcomes`, even when
they have no reconstruction to replay. Missing original inputs cause import
to fail, rather than reacquiring or treating the event as zero support.

```python
from isovar.audit_bundle import AuditBundle
from isovar.sid_audit_bundle import bundle_sid_audit
from isovar.audit_replay import replay_sv_rna

bundle = bundle_sid_audit("EXISTING_AUDIT", "saved-sid-audit")
# After moving or restoring the directory, reopen its manifest.
bundle = AuditBundle.open("relocated-sid-audit")
run_id = next(iter(bundle.manifest["metadata"]["sv_rna_runs"]))
checkpoint = replay_sv_rna(bundle, run_id, "replay-checkpoints")
assert checkpoint["comparison"] == "matched"
```

Import preserves inventory, transcript-reference snapshots, receipts, source
headers, regional BAMs/indexes, shared provenance assets and reconstruction
checkpoints. Original file paths remain provenance; runtime paths come from
logical bundle names. Frozen transcript sequences/exons/CDS are enough for the
existing SV reconstruction engine, so replay does not require a local Ensembl
installation. Add annotation files required by other analyses explicitly with
`extra_files={name: path}` or repeated `--extra-file NAME=PATH`. Other analysis
files can be bundled and verified through the generic API; this replay command
implements SV RNA reconstruction, not every downstream report or small-variant
analysis.

Replay verifies the whole bundle, uses recorded read filters and reconstruction
parameters, preserves acquisition limitations, and compares the **entire**
result with the frozen reconstruction. It writes a new checkpoint outside the
bundle. Repeating the command reuses only a matching, intact checkpoint.
Software changes fail by default; `allow_software_change=True` or
`--allow-software-change` permits an explicit new-engine comparison. A differing
result still fails and preserves its checkpoint for inspection. Imported old
software identities and wall-time budgets remain in `original_request`; replay
uses the bundled settings with the recorded bundle-creation software, without
the original runner's process-level alarm. Use your job runner for wall-time
limits. `scratch_dir` / `--scratch-dir` only selects disposable dense-engine
scratch storage. Sequence, support and biological interpretation are unchanged.

## Manifest and errors

For other datasets, `AuditBundle.create` accepts replay plans in
`metadata["sv_rna_runs"]`, keyed by your own stable run IDs. Each plan requires
`bam`, `index`, `input` (logical names of bundled files), `source` (the original
alignment identity), and `parameters` (keyword arguments of
`reconstruct_sv_rna`). The input JSON uses the same public schema as
[`isovar sv-rna`](sv-rna.md). Optional `read_collection` contains `ReadCollector`
keyword arguments. Optional `expected` names a frozen full reconstruction JSON;
without it, comparison is explicitly `not_requested`. `input_limitations`,
`upstream_acquisition` and `original_request` preserve acquisition/provenance.
Misspelled or missing plan fields are errors. A minimal plan is:

```json
{
  "sv_rna_runs": {
    "tumor-T2-event-forward": {
      "bam": "reads/tumor.bam",
      "index": "reads/tumor.bam.bai",
      "input": "events/event.json",
      "source": "original-library-T2",
      "parameters": {"assemble": false, "dense_support": true},
      "expected": "results/original.json"
    }
  }
}
```

Pass this object as `metadata` and declare each named file in `files`.

`manifest.json` has `format="isovar-audit-bundle"`, `schema_version=1`,
`files`, `metadata`, `software` and `bundle_id`. Each file is stored at
`files/LOGICAL_NAME`; its manifest entry has `sha256` and `size`. `bundle_id` is
SHA-256 of canonical sorted JSON excluding itself. Unknown schema versions,
escaping paths, missing files, changed bytes and incomplete archives fail
validation. An embedded manifest verifies consistency; it is not an external
signature. Compare `bundle_id` with your trusted recorded identity when
receiving a bundle.

`AuditBundle.open()` verifies every file by default. `verify()` returns the
identity, verified file count and payload bytes. `manifest` returns a detached
JSON dictionary. The supported APIs are `AuditBundle`, `AuditBundleError`,
`bundle_sid_audit` and `replay_sv_rna` in their documented modules. Validation
errors raise `AuditBundleError`; ordinary filesystem errors retain their Python
exception types. The CLI prints JSON on stdout; errors go to stderr with exit
status 2. It never executes commands supplied by a bundle.
