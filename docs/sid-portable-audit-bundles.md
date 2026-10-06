# Preserved Sid audit bundles

The installed [audit-bundle API and CLI](audit-bundles.md) preserve the selected
original inputs, frozen references, receipts, screen/report artifacts and SV
reconstruction checkpoints in three independently movable bundles. Original
receipt/result bytes are unchanged; runtime files are resolved from logical
names in the bundle. The compressed manifests below are checked into Git.
The regional reads and annotation archives remain in durable local storage.

| Bundle | Verified files | Payload bytes | SV replay views | Frozen manifest |
| --- | ---: | ---: | ---: | --- |
| Vaccine/control selection | 239 | 556,079,597 | 40 | [Manifest](data/sid-portable-bundles/priority.manifest.json.gz) |
| Non-SNV/fusion additions | 535 | 3,858,891,048 | 16 | [Manifest](data/sid-portable-bundles/additions.manifest.json.gz) |
| Corrected dense controls | 34 | 1,376,437,179 | 4 | [Manifest](data/sid-portable-bundles/dense.manifest.json.gz) |

The additions include all eight pinned Ensembl annotation files used by the
124 allele/product screens. Repeated references to the same source file are
validated once during import; copying still checks the expected source checksum
and the copied bytes. No whole RNA library was downloaded for this export.

The durable directories/archives are under this checkout's `.cache/`:

```text
portable-sid-priority-final-20261006/
portable-sid-priority-final-20261006.tar.gz
portable-sid-additions-20261006/
portable-sid-additions-20261006.tar.gz
portable-sid-dense-final-20261006/
portable-sid-dense-final-20261006.tar.gz
```

Copy an archive to your backup storage using ordinary tools. After relocation,
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
full-ORF supporting fragments. All 5,794 tests and lint pass. A stale module in
an existing developer `build/lib` directory was detected by software identity
checks ([#454](https://github.com/openvax/isovar/issues/454)); final wheel
validation uses a clean checkout, matching the canonical deployment workflow.

All 60 completed SV views retain their original acquisition/outcome provenance.
Other unfinished or unselected catalogue targets remain outside this saved
scope. The importer preserves small-variant and report assets, while the public
replay command implements SV RNA reconstruction; downstream report/small-variant
commands retain their separately documented input/annotation requirements.
