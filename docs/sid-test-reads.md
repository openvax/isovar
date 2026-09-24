# Minimal Sid test reads

Isovar's tests use original reads from the public Sid osteosarcoma data. Only
the records the tests need ship with the package, so the tests run offline.
This page covers listing and exporting those reads, regenerating them from the
source BAMs, and changing which reads are selected.

Two test-data sets exist; this page is about the first:

| Data | Used by | Contents |
|---|---|---|
| This bundle | Isovar's own tests | The exact SAM records each test needs, shipped in the wheel |
| [Shared osteosarc data](osteosarc-data.md) | Vaxrank, and anyone needing BAM files | A pinned 49-case set of BAM and index files, downloaded to a shared cache |

## What is in the bundle

- `sid-reads.json.gz` stores each record once. `sid-read-recipe.json.gz` lists
  311 ordered fixture selections from 39 source assets, the tests that use them,
  and the bounded intervals they were taken from.
- There are 27,077 unique source/SAM-record pairs. Reusing a record in several
  tests costs no extra storage, and genuine duplicate records are restored on
  export.
- The selections cover SNVs and indels, native GRCh37 and mitochondrial loci,
  long reads, fusions, phasing, clipping and matched DNA/RNA regressions. The
  external K562 control is kept separately.
- They are test selections, not a random sample or a whole locus, and cannot
  estimate sample coverage or VAF. Historical allele definitions and protein
  expectations are unchanged.

Acquisition, record selection, transport and integrity checks are handled by the
installed `osteosarc` package (0.2.3 or later), through snapshot Assets and
indexed Region extraction. The shipped bundle was generated with osteosarc 0.1.0,
which its recipe records; its read bytes are unchanged since. New recipe
compilations record the installed osteosarc version.

## Use the installed package offline

No network access is needed to list, verify or export:

```sh
python -m isovar.sid_data verify
python -m isovar.sid_data list
python -m isovar.sid_data export \
  --fixture 'fusions/long-read/FOXO3--STRADA-CCDC47.ONT-T1-tagged.sam.gz' \
  --output /tmp/foxo3-test-reads.sam
```

In Python, `isovar.sid_data` has `verify()` and `export_fixture(name, output)`.

Export keeps every SAM field and tag, the qualities (including absent QUAL), the
original record order and repeated records. It refuses to overwrite an existing
file. The reduced header keeps the reference identities, the read groups the
reads refer to, and the program records and predecessor chains those reads
need. It is marked `SO:unknown`, because some selections order reads
deliberately rather than by coordinate. The exported SAM matches the original
record for record, but is not byte-identical to the historical BAM files and
their full headers.

## Regenerate from osteosarc

Regeneration is explicit and needs samtools built with HTTPS support. A new
cache needs a metadata snapshot first:

```sh
python -m osteosarc --cache /tmp/sid-cache sync isovar-test-reads
python -m isovar.sid_data generate \
  --snapshot isovar-test-reads --cache /tmp/sid-cache \
  --acquired /tmp/sid-selected
python -m isovar.sid_data pack \
  --acquired /tmp/sid-selected --output /tmp/sid-reads.json.gz
python -m isovar.sid_data verify --bundle /tmp/sid-reads.json.gz
```

How generation checks its inputs:

- Each source Asset's ID, URL, bucket key, size and modification time must match
  the recipe. A changed asset, or a missing or modified required record, is an
  error that needs review. An existing named snapshot can be reused, and a newer
  one is acceptable only when the pinned asset identities still match.
- The recipe records the original snapshot ID. Each acquisition also records the
  snapshot actually used, the remote identity, index and header receipts,
  resolved regions and hashes.
- Queries use one reference position per required alignment, merging adjacent
  positions. osteosarc returns whole original records overlapping those
  intervals, and the generator keeps **only the recipe's exact SAM checksums**,
  including the required number of identical source records. References to one
  alignment from several fixtures are told apart from true duplicates.
- No sequence is cropped or realigned. Mates and supplementary records a test
  needs have their own query positions. Adding unmapped reads requires a
  reviewed mate-recovery recipe, not a whole-file scan.

Coordinates follow the [SAM specification](https://samtools.github.io/hts-specs/SAMv1.pdf)
and [samtools region semantics](https://www.htslib.org/doc/samtools-view.html).

Other generation options:

- The acquisition directory is resumable: completed selections are verified
  before reuse.
- `--source ASSET_ID_PREFIX` limits a run to one source in the recipe.
- `--offline` forbids network access and requires cached regional extracts.
  The cache can hold larger temporary regional BAMs and indexes; these never go
  into the package, and no test triggers acquisition.
- For a shared, versioned osteosarc panel recipe, `generate` also accepts
  `--panel-recipe recipe.json --panel-source SOURCE_ID=original.bam --output
  NEW_DIRECTORY --offline`, with one `--panel-source` per local original input.
  Panel output contains indexed BAMs, checksums, record multiplicities, source
  and header identities, and selection reasons. See osteosarc's
  [shared fixture workflow](https://github.com/iskandr/osteosarc/blob/v0.2.3/docs/fixture-migration.md)
  and [Osteosarc #15](https://github.com/iskandr/osteosarc/issues/15).

Write each regeneration to a new destination. The checked-in scientific
expectations remain the reference.

## Change the selected reads

In a repository checkout, the reviewed selection inventory is
[`tests/data/osteosarc/bundle.py`](../tests/data/osteosarc/bundle.py). Its
entries connect source fixtures to the tests that use them, including original
SAM records embedded in small JSON evidence fixtures. Compile a new recipe with:

```sh
python -m tests.data.osteosarc.bundle \
  --snapshot isovar-test-reads --cache /tmp/sid-cache \
  --output /tmp/sid-read-recipe.json.gz
```

Then:

1. pass `--recipe /tmp/sid-read-recipe.json.gz` to `generate`, `pack` and `verify`;
2. review the selection diff and regenerate through osteosarc;
3. run `tests/test_sid_data.py` and the tests that use the changed fixtures
   before replacing the packaged data.

The checks require every packaged record to belong to a listed selection, and
they compare complete SAM strings, not just names or counts.

The regional audit builders also use osteosarc. Configure them with
`ISOVAR_SID_SNAPSHOT` and, optionally, `ISOVAR_SID_CACHE`; they can produce
larger research inputs outside the package. The package excludes those inputs,
full source BAMs, index caches, transcript and reference archives, audit and
benchmark reports, and figures. Scientific report and reference fixtures stay in
the repository. The original
container checksums remain historical regression checks, and the packaged
bundle is checked against the same original records.

The read data are CC0-1.0. Source attribution and licensing are kept in the
recipe and in [`source_registry.yaml`](../tests/data/osteosarc/source_registry.yaml).
