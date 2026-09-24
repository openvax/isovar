# Shared osteosarc regression data

Isovar ships the manifest of a fixed set of real BAM files from the public Sid
osteosarcoma data, `osteosarc / vaccine-rna-v1`, which Vaxrank also uses for its
regression tests. The files themselves are downloaded to a shared cache on
request. For the individual SAM records that Isovar's own tests use, which ship
inside the package, see [Minimal Sid test reads](sid-test-reads.md).

## What the manifest pins

- 49 cases covering 44 original vaccine loci, with 98 BAM and index files
  totaling 4,126,730 bytes.
- Each file's size and SHA-256 checksum, and each case's native alleles, source
  product and selection provenance. The manifest also pins the Isovar commit that produced
  it, `0cad5b275c852263a1c77722aa463fe5568b2c76`.
- Native GRCh37 NR2F2 and mitochondrial cases, with their original reference
  identities.

These are deliberately selected regression reads, not a regional analysis or a
coverage estimate. Correcting catalogue entries and choosing fixtures are
separate workflows. This revision keeps its original MAP2 allele: adopting
osteosarc's corrected complex allele needs a separate biological review
([Vaxrank #489](https://github.com/openvax/vaxrank/issues/489)). The pinned
manifest and BAM files do not change when osteosarc updates its catalogue.

## Import, download and export

Acquisition goes through the cache API of the installed `osteosarc` package:

```sh
# Import the fixtures already in this repository into the shared cache, offline.
python -m isovar.osteosarc_data --offline \
  --import-corpus tests/data/osteosarc/expansion/corpus

# Or download the pinned files into the shared cache.
python -m isovar.osteosarc_data

# Write a new offline copy from the cache.
python -m isovar.osteosarc_data --offline --output /path/to/new/fixture-export

# Verify an export without touching the cache.
python -m isovar.osteosarc_data --verify-only --output /path/to/fixture-export

# Refetch corrupt cache files; valid ones are reused.
python -m isovar.osteosarc_data --repair-cache
```

An export contains the portable manifest and the exact original BAM and index
bytes. It performs no alignment, filtering or further sampling. Ordinary tests
read their checked-in fixtures and never download anything.

The fuller scientific selection records and pinned transcript references are in
`tests/data/osteosarc/expansion/corpus`, and its `fixtures.py` recipe regenerates
them from regional evidence. Run it into a new directory and compare the results
before proposing a new dataset revision. Exporting this revision does not require
an Isovar version bump.

## The cache

The default is osteosarc's shared OpenVax cache. Choose another root with
`--cache-root`, `OSTEOSARC_CACHE` or `OPENVAX_DATA_CACHE`, checked in that order;
otherwise the platform's cache directory for `openvax` is used. Files are stored
as:

```text
<root>/objects/sha256/<sha256><original suffixes>
```

This is the same layout as Vaxrank's downloader. Valid files already there are
verified and adopted, with osteosarc import receipts, without downloading or
rewriting them. Receipt timestamps record the local import time, and the
manifest keeps the original provenance. Consumers need no Vaxrank dependency:
the manifest and the osteosarc cache are enough.

## Failures

- Every file must match both its size and its checksum. Missing files offline,
  and corrupt files, fail visibly.
- Replacing corrupt cache files needs `--repair-cache`, which goes online and
  cannot be combined with `--offline`.
- An existing export is verified, never repaired or overwritten, even when it is
  empty or partial. Every input is resolved before a new export starts, and its
  manifest is written last, so an interrupted export fails verification. Retry
  into a new directory.

The implementation is [`isovar/osteosarc_data.py`](../isovar/osteosarc_data.py),
and the packaged manifest is
[`isovar/data/osteosarc-test-data-v1.json`](../isovar/data/osteosarc-test-data-v1.json).
Its tests check complete SAM-record multisets (including qualities, flags and
tags), native case provenance, reuse of the shared cache, offline failures and
explicit repair.
