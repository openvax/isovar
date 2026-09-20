# Minimal Sid test reads

Isovar's wheel and source distribution bundle only the original Sid read records
selected by its tests. `sid-reads.json.gz` stores shared records once;
`sid-read-recipe.json.gz` specifies 311 ordered fixture selections from 39 source
assets, their consuming tests, and bounded acquisition intervals. There are
27,077 unique source/SAM-record pairs. Reusing a record across tests does not
increase storage, while actual duplicate multiplicity is preserved on export.

The bundle is generated through **osteosarc 0.1.0**, using snapshot Assets and
indexed Region extraction. It includes SNV/indel, native GRCh37 and mitochondrial,
long-read, fusion, phasing, clipping, and matched DNA/RNA regression inputs.
The external K562 control remains separate. These are test selections, not a
random sample, a whole locus dataset, or an estimator of sample coverage/VAF.
Historical allele definitions and protein expectations are unchanged.

## Use the installed package offline

No osteosarc installation or network access is needed to list, verify or export:

```sh
python -m isovar.sid_data verify
python -m isovar.sid_data list
python -m isovar.sid_data export \
  --fixture 'fusions/long-read/FOXO3--STRADA-CCDC47.ONT-T1-tagged.sam.gz' \
  --output /tmp/foxo3-test-reads.sam
```

The same operations are available as `verify()` and
`export_fixture(name, output)` in `isovar.sid_data`. Export preserves every
SAM field, tag, quality (including absent QUAL), original record order, and
repeated records. It refuses to overwrite an existing file. The reduced header
retains reference identities, referenced read groups, and the program records
and predecessor chains needed by those reads. It is marked `SO:unknown` because
some test selections deliberately order observations rather than coordinates.
Exported SAM is equivalent at the record level, not byte-identical to the
historical BAM containers and their full headers.

## Regenerate from osteosarc

Regeneration is explicit and needs Python 3.10+, `isovar[data]`, and HTTPS-enabled
samtools. A new cache needs a metadata snapshot first:

```sh
python -m pip install 'isovar[data]'
python -m osteosarc --cache /tmp/sid-cache sync isovar-test-reads
python -m isovar.sid_data generate \
  --snapshot isovar-test-reads --cache /tmp/sid-cache \
  --acquired /tmp/sid-selected
python -m isovar.sid_data pack \
  --acquired /tmp/sid-selected --output /tmp/sid-reads.json.gz
python -m isovar.sid_data verify --bundle /tmp/sid-reads.json.gz
```

An existing named snapshot may be reused. Generation compares each exact source
Asset's ID, URL, bucket key, size and modification time to the recipe. A changed
asset or missing/modified required record is an error requiring review. The
recipe records the original snapshot ID; every acquisition also records the
snapshot actually used, remote identity, index/header receipts, resolved
regions, and hashes. A new metadata snapshot is acceptable only when the pinned
source asset identities still match.

Queries use one reference position per required alignment, coalescing adjacent
positions. osteosarc retrieves whole original records overlapping those indexed
intervals; the generator then retains **only the recipe's exact SAM checksums**.
It verifies the required number of identical source records as well. Fixture
cross-references to the same alignment are distinguished from true source
duplicates. No sequence is cropped or realigned; mates and supplementary records
needed by a test have their own explicit query positions. Unmapped additions
require a reviewed mate-recovery recipe rather than a whole-file scan.
Coordinates follow the [SAM specification](https://samtools.github.io/hts-specs/SAMv1.pdf)
and [samtools region semantics](https://www.htslib.org/doc/samtools-view.html).

The acquisition directory is resumable: completed selections are verified
before reuse. `--source ASSET_ID_PREFIX` restricts a generation run to a source
listed in the recipe. `--offline` forbids network and requires cached regional
extracts. A cache may contain larger temporary regional BAMs and indexes; those
are never copied into the package. No test triggers acquisition.

## Change the selected tests deliberately

From a repository checkout, the reviewed selection inventory is
[`tests/data/osteosarc/bundle.py`](../tests/data/osteosarc/bundle.py). Its explicit
entries connect source fixtures to consuming tests, including original SAM
records embedded in small JSON evidence fixtures. Compile a new recipe with:

```sh
python -m tests.data.osteosarc.bundle \
  --snapshot isovar-test-reads --cache /tmp/sid-cache \
  --output /tmp/sid-read-recipe.json.gz
```

Pass `--recipe /tmp/sid-read-recipe.json.gz` to `generate`, `pack` and `verify`.
Review the selection diff, regenerate through osteosarc, and run
`tests/test_sid_data.py` plus the consuming tests before replacing packaged data.
The inventory/checks require every packaged record to belong to a listed test
selection; they compare complete SAM strings, not merely names or counts.

The existing regional audit builders also use osteosarc. Configure them with
`ISOVAR_SID_SNAPSHOT` and optionally `ISOVAR_SID_CACHE`. They can still produce
larger research audit inputs outside the package. The package excludes those
inputs, full source BAMs, index caches, transcript/reference archives, audit and
benchmark reports, and figures. Scientific report/reference fixtures remain in
the checkout. Their original container checksums remain historical regression
checks; the packaged read bundle is checked against the same original records.

The read data are CC0-1.0; source attribution and licensing are retained in the
recipe and [`source_registry.yaml`](../tests/data/osteosarc/source_registry.yaml).
The earlier [shared BAM cache workflow](osteosarc-data.md) remains available for
consumers needing the 49-case `vaccine-rna-v1` manifest and exact BAM/index bytes.
