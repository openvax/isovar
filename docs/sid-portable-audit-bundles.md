# Preserved Sid audit bundles

The installed [audit-bundle API and CLI](audit-bundles.md) preserve the selected
original inputs, frozen references, receipts, screen/report artifacts and SV
reconstruction checkpoints in independently movable bundles. Original
receipt/result bytes are unchanged; runtime files are resolved from logical
names in the bundle. The compressed manifests below are checked into Git.
The regional reads and annotation archives remain in durable local storage.
The three historical backups below total **2,172,408,842 bytes (2.02 GiB)**; exact archive
sizes, SHA-256 and bundle identities are in [the archive index](data/sid-portable-bundles/archives.json).

| Bundle | Verified files | Payload bytes | SV replay views | Frozen manifest |
| --- | ---: | ---: | ---: | --- |
| Vaccine/control selection | 239 | 556,079,597 | 40 | [Manifest](data/sid-portable-bundles/priority.manifest.json.gz) |
| Non-SNV/fusion additions | 535 | 3,858,891,048 | 16 | [Manifest](data/sid-portable-bundles/additions.manifest.json.gz) |
| Corrected dense controls | 34 | 1,376,437,179 | 4 | [Manifest](data/sid-portable-bundles/dense.manifest.json.gz) |

The additions include all eight pinned Ensembl annotation files used by the
124 allele/product screens. Repeated references to the same source file are
validated once during import; copying still checks the expected source checksum
and the copied bytes. No whole RNA library was downloaded for this export.

The durable archives are under this checkout's `.cache/`. Regenerable export
and validation-restore directories were reclaimed after every payload byte was
verified in a retained archive; original audit inputs and manifests remain:

```text
portable-sid-priority-final-20261006.tar.gz
portable-sid-additions-20261006.tar.gz
portable-sid-dense-final-20261006.tar.gz
```

Copy an archive to your backup storage using ordinary tools. The three historical
bundles pin Isovar 1.45.2; use that version to reproduce their original software
identity. After relocation,
restore it to a new directory, compare its identity with the frozen manifest,
inspect its run IDs and replay a selected SV view:

```sh
isovar audit-bundle unpack portable-sid-priority-final-20261006.tar.gz \
  --output RESTORED
isovar audit-bundle verify RESTORED
isovar audit-bundle inspect RESTORED
isovar audit-bundle replay RESTORED --run-id SOURCE/EVENT/reverse --output REPLAY
```

Installed-wheel validation ran from `/private/tmp` with Python isolated mode,
using only the built wheel and ordinary dependencies. It restored the archive,
verified every file, reproduced the entire frozen ONT T2 TPST1–CRCP reverse-view
reconstruction with Python networking blocked, and reused the matching
checkpoint through the installed CLI. Its 30-aa candidate
`MPSRRRGGSSLIPGRMPILRFSVTTRYFSY` retains **172 full-interval fragments / 139
barcode labels**. Counts retain their original completeness/lineage flags;
this verifies RNA sequence/support, not translation. The exact result hash and
validation are in [the replay ledger](data/sid-portable-bundles/wheel-replay.json).

The offline regression also restores real original PARD3B SAM reads after
deleting its original directory and verifies the expected sequence and four
full-ORF supporting fragments. A stale module in
an existing developer `build/lib` directory was detected by software identity
checks ([#454](https://github.com/openvax/isovar/issues/454)); final wheel
validation uses a clean checkout, matching the canonical deployment workflow.

All 60 completed SV views retain their original acquisition/outcome provenance.
Other unfinished or unselected catalogue targets remain outside this saved
scope. The importer preserves small-variant and report assets, while the public
replay command implements SV RNA reconstruction; downstream report/small-variant
commands retain their separately documented input/annotation requirements.

## Dense recovery publication

The fourth archive, `portable-sid-recovery-final-20261011.tar.gz`, preserves
all 44 selected recovery outcomes, with 42 executable replay plans. The ONT
IGKC timeout and interrupted scan retain unknown support; they do not have
replay plans. [Follow-up #470](https://github.com/openvax/isovar/issues/470)
tracks their separately produced successors below.

The [producer-wheel proof](data/sid-portable-bundles/recovery-wheel-replay.json)
verifies every archived payload byte, relocation of the full bundle, installed
CLI restoration of a complete selected view, offline reconstruction and CLI
checkpoint reuse. Both wheel checks matched the same complete result hash.

It contains **722 verified files / 11,348,954,838 payload bytes**.
The archive is **7,043,781,669 bytes**; its SHA-256 and bundle identity
are in the [archive index](data/sid-portable-bundles/archives.json) and
[frozen manifest](data/sid-portable-bundles/recovery.manifest.json.gz).
The [report dependency ledger](data/sid-portable-bundles/recovery-report-dependencies.json.gz)
accounts for all 284 frozen inputs across the four declared bundles. Preserve
all four archives to retain the combined report's original inputs.

The recovery producer was an **unreleased 1.45.3 development build**. Its
preserved wheel is `.cache/sid-recovery-wheel-20261006/isovar-1.45.3-py3-none-any.whl`;
it is not a PyPI release. The installed 1.47.1 wheel also reproduced the complete
selected ONT CDK11 forward view after archive restoration with Python networking
blocked. The comparison explicitly allowed the software change and matched the
entire frozen result; see the [release-wheel proof](data/sid-portable-bundles/recovery-release-wheel-replay.json).
With 1.47.1, restore and inspect the recovery archive using the commands above,
then replay a selected run with `--allow-software-change`. A changed result
still fails comparison.

## ONT IGKC completion

Two new data archives preserve the independently verified forward and reverse
ONT IGKC successors, their native and normalized reconstruction bytes, original
regional BAM/index, frozen references, producing wheel, previous failures and
validation controls. A third addition contains the refreshed reports and their
provenance. The preceding four archives and bundle identities are unchanged.

- `portable-sid-ont-forward-completion-20261014.tar.gz`: 198 files, 21,239,469,931 payload bytes; 21,222,114,870 archive bytes. [Manifest](data/sid-portable-bundles/ont-forward-completion.manifest.json.gz).
- `portable-sid-ont-reverse-completion-20261014.tar.gz`: 237 files, 14,193,369,864 payload bytes; 14,165,616,980 archive bytes. [Manifest](data/sid-portable-bundles/ont-reverse-completion.manifest.json.gz).
- `portable-sid-ont-completion-reports-20261014.tar.gz`: 128 files, 393,109,493 payload bytes; 386,332,536 archive bytes. [Manifest](data/sid-portable-bundles/ont-completion-reports.manifest.json.gz).

The [archive index](data/sid-portable-bundles/archives.json) pins all seven archives.
The [dependency inventory](data/sid-portable-bundles/ont-completion-report-dependencies.json)
resolves every input of the refreshed screen and event report to its exact file
in those archives. Preserve all seven for the combined report. Each orientation
archive contains everything needed to replay its own view.

Both data archives and the report addition were actually restored and verified.
The [installed-wheel proof](data/sid-portable-bundles/ont-completion-installed-replay.json)
records full forward reconstruction from a relocated restore with Python
networking blocked, matching every frozen result value and scanning all
2,715,459 input records. CLI resume verifies the complete compressed checkpoint
and leaves its bytes unchanged. The reverse archive was restored and inspected;
its producing result passed independent sequence, support and translation checks.

The replay wheel includes the streamed decoder and checkpoint validator; its
identity is recorded separately from the frozen producing wheel. Duplicate
working files and the disposable replay checkpoint were reclaimed only after
verification, with checksummed cleanup receipts. Complete native and normalized
scientific evidence remains in the archives. These checks preserve the recorded
acquisition/discovery limitations and do not establish translation, tumor
specificity or an independent DNA rearrangement.
