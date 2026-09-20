# Shared osteosarc regression data

For the packaged test-only read selections and regeneration from source BAMs,
see [Minimal Sid test reads](sid-test-reads.md). This page covers the older
shared 49-case BAM/index cache, which remains available unchanged.

Isovar packages the `osteosarc / vaccine-rna-v1` manifest used by Vaxrank:
49 cases covering 44 original vaccine loci, with 98 BAM/index objects totaling
4,126,730 bytes. The manifest pins Isovar commit
`0cad5b275c852263a1c77722aa463fe5568b2c76`, file sizes, SHA-256 checksums,
native alleles, source products and selection provenance. It includes native
GRCh37 NR2F2 and mitochondrial cases with their original reference identities.

These are deliberately selected regression reads. Full regional analysis,
source coverage estimates, live catalogue corrections and fixture selection
remain separate workflows. In particular, this historical revision retains
its MAP2 allele; adopting osteosarc's corrected complex allele requires a
separate biological review ([Vaxrank #489](https://github.com/openvax/vaxrank/issues/489)).
The pinned manifest and BAMs do not change when osteosarc updates its catalogue.

## Import, download and export

Acquisition uses the published `osteosarc==0.1.0` cache API. It requires
Python 3.10+ and the optional data extra:

```sh
python -m pip install 'isovar[data]'

# Import the existing Isovar fixtures into the shared cache without network.
python -m isovar.osteosarc_data --offline \
  --import-corpus tests/data/osteosarc/expansion/corpus

# Or explicitly download the pinned objects into the shared cache.
python -m isovar.osteosarc_data

# Materialize a new offline copy from the cache.
python -m isovar.osteosarc_data --offline --output /path/to/new/fixture-export

# Verify an export without touching the cache or requiring osteosarc.
python -m isovar.osteosarc_data --verify-only --output /path/to/fixture-export

# Explicitly refetch corrupt cache objects; valid objects are reused.
python -m isovar.osteosarc_data --repair-cache
```

Isovar's core Python 3.9 support is unchanged. Importing the module, displaying
help and verifying existing exports do not require osteosarc. Ordinary tests
read their checked-in fixtures and perform no acquisition.

An export contains the portable manifest and exact original BAM/index bytes.
It performs no alignment, filtering or further sampling. The fuller scientific
selection records and pinned transcript references remain in
`tests/data/osteosarc/expansion/corpus`; the existing `fixtures.py` recipe owns
regeneration from regional evidence. Run it into a new directory and compare
the results before proposing a new dataset revision. Exporting this revision
is independent of an Isovar package-version bump.

## Cache reuse and failures

The default is osteosarc's shared OpenVax cache. Select a root with
`--cache-root`, `OSTEOSARC_CACHE`, or `OPENVAX_DATA_CACHE`, in that order.
Otherwise the platform cache directory for `openvax` is used. Objects live at:

```text
<root>/objects/sha256/<sha256><original suffixes>
```

This is the same layout as Vaxrank's downloader. Valid objects already there
are verified and adopted with osteosarc import receipts without redownloading
or rewriting them. Receipt timestamps record local import time; the manifest
retains original source/selection provenance. No Vaxrank Python dependency is
needed. Consumers can use the manifest and osteosarc cache directly.

Every object must match both size and checksum. Missing offline objects and
corrupt objects fail visibly. Online `--repair-cache` is required to replace
corrupt cache bytes; it cannot be combined with `--offline`. Existing exports
are verified and never repaired or overwritten, including empty or partial
directories. All inputs are resolved before creating a new export, and its
manifest is written last. An interrupted export is rejected on verification;
choose a new output directory to retry.

The implementation is in [`isovar/osteosarc_data.py`](../isovar/osteosarc_data.py),
and the packaged manifest is
[`isovar/data/osteosarc-test-data-v1.json`](../isovar/data/osteosarc-test-data-v1.json).
The regression suite verifies complete SAM-record multisets (including
qualities, flags and tags), native case provenance, shared-cache reuse,
offline failures and explicit repair.
