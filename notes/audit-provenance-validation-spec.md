# Audit provenance validation

Validate the evidence recorded by existing Sid reports before producing a new
portable manifest. Named `pins` maps and `result-pins.json.gz` map audit-relative
paths to SHA-256 values; structured file pins retain their existing document
base. Reject missing, changed or malformed map entries without publishing a
bundle. Check acquisition, source and root identities for every checkpoint,
including outcomes without a reconstruction; only replay generation is optional.

Use one 16 MiB manifest limit for creation, opening, verification, packing and
restoration. Count encoded bytes, including the writer's final newline. Failed
creation or backup must leave no published destination. Exercise the exact
boundary, oversized metadata, on-disk whitespace and archive restoration.

Add regression tests for both pin-map formats, modified checkpoints, unresolved
outcomes and backup/restore consistency. Validate with the latest mutually
compatible stable release dependencies from live PyPI, record the resolved
versions and artifact hashes, and run lint, the full suite and GitHub CI. Keep
historical producer/reference identities immutable. Ship via a versioned PR and
the standard deployment script.

## Release dependencies

The [live dependency check](data/audit-validation-dependencies.json) records the
Python 3.14.6 macOS ARM64 release environment, including the downloaded artifact
SHA-256 values checked against PyPI. All 70 dependencies and release tools matched
their latest non-yanked stable releases on 2026-10-14. This is a validation
snapshot, not a cross-platform installation lock. That validation's CI tested
Python 3.9, 3.10 and 3.11 with the newest compatible releases available to each
interpreter; current Isovar releases require Python 3.11 or later.

The release environment used PyEnsembl 2.27.3, Varcode 10.15.0,
Orfdata 0.18.0 and osteosarc 0.15.5. Existing version ranges retained older
Python support at that time. The completed ONT producer and installed replay
retain their exact original environments. Their [frozen replay dependency
catalogue](data/sid-sv-recovery/ont-completion/frozen-replay-dependencies.json)
is preserved separately from these current release pins. Historical scientific
evidence pins continue to identify their original bytes.

GitHub Actions' live releases were checked as well: checkout v7.0.1,
setup-python v7.0.0 and cache v6.1.0 are covered by the workflow's existing
v7/v7/v6 major tags.
