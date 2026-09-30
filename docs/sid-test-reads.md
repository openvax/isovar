# Sid test reads

Isovar's tests use original reads from the public Sid osteosarcoma data (CC0).
They ship with Isovar, in `isovar/data/sid-reads`, so the tests never need
network access. That folder is an osteosarc fixture bundle holding Isovar's
members of **openvax-v1**, the test-data bundle osteosarc publishes for all four
OpenVax libraries (Isovar, Varcode, Vaxrank and Topiary). osteosarc makes it
from openvax-v1, and CI checks that it still does, record for record.

## What Isovar's bundle holds

- **345 legacy members, named `isovar/<path>`** after the files Isovar's tests read
  (`<path>` is relative to `tests/data`). They hold exactly the same records.
- **124 whole read files**, 5.8 MB of SAM and BAM:
  - the six-locus RNA;
  - the 49-case vaccine corpus and its primary-only copies;
  - the stress corpus;
  - the PIP5K1A and PacBio comparisons;
  - the long-read fusions;
  - the chimeric examples.
- **217 selections referenced by JSON fixtures**: `file#/pointer` names the
  original location of the records. The 25 JSON files now store member names
  and SAM-text digests instead of duplicate SAM lines.
- **4 single-read selections**, named `read_ends/osteosarc.sam#READ`, consumed
  directly from the bundle. The former SAM copy is removed.
- **184 T2 short-read members** retain their upstream names and original
  selection, bringing the packaged total to **529 members**. They cover every
  literal GRCh38 small-variant target in openvax-v1. The complete inventory also
  retains three GRCh37/hg19 targets without a corresponding T2 member and two
  unresolved catalog alleles. The [full-panel regression](../tests/data/osteosarc/union/README.md)
  accounts for all of them and pins independent allele and protein checks.
- **Selection and provenance.** osteosarc owns read selection and records the
  source, snapshot and reason for every template. openvax-v1 also covers all
  187 osteosarc site variants and the other libraries' fixtures. See
  [osteosarc's shared test data](https://iskandr.github.io/osteosarc/test-data/#shared-test-data-openvax-v1)
  and [iskandr/osteosarc#56](https://github.com/iskandr/osteosarc/issues/56).

These are test selections, not a random sample or a whole locus, and they
can't estimate sample coverage or VAF. Historical allele definitions and
protein expectations are unchanged.

## Using the reads

`isovar.sid_data` checks the packaged bundle against its pinned manifest hash
(`PACKAGED_MANIFEST_SHA256`). osteosarc's `bundle_file` exports each read file
from it once, in its original format, and later runs reuse the export:

```python
from isovar import sid_data

path = sid_data.path("osteosarc/bulk_star_t0.sam.gz")   # .../bulk_star_t0.sam.gz.sam.gz
sid_data.export("fusions/corpus/TPST1--CRCP-T1.input.json.gz#/original_records", "reads.sam")
sid_data.members()                                      # every Isovar member
records = sid_data.records("chimeric/osteosarc-ont.sam")  # SAM digest -> SAM text
fixture = sid_data.read_fixture("tests/data/fusions/corpus/TPST1--CRCP-T1.input.json.gz")
```

```sh
python -m isovar.sid_data list
python -m isovar.sid_data path osteosarc/expansion/corpus/28-NTF3-chr12-5494381-53f498a544883d51.bam
python -m isovar.sid_data export chimeric/osteosarc-ont.sam --output /tmp/chimeric.sam
```

- **Size.** The bundle is about 37 MB, mostly JSON (its manifest, recipe and the
  sources' original headers).
- **Where the exports live.** In osteosarc's cache: `OSTEOSARC_CACHE`, or the
  shared OpenVax cache (`OPENVAX_DATA_CACHE`, or the platform's `openvax`
  cache directory). osteosarc exports every member of each source BAM it
  reads, so the exports take about 35 MB.
- **Names.** Exported files are read-only and named `<member>.<format>`, such
  as `bulk_star_t0.sam.gz.sam.gz`; `path` returns them.
- **A missing, changed or damaged bundle** (even one with an extra file, such
  as a `.DS_Store`) raises `SidDataUnavailable`, which says to reinstall Isovar
  or build the bundle again.

An exported file holds exactly the member's original records, but it is
coordinate-sorted and carries the source's full header. So its bytes, and the
order of records at the same position, differ from the old copies. Tests
compare records, not file checksums. `tests/test_sid_data.py` checks every
decoded selection against the bundle with osteosarc's `check_fixtures`.
It also checks fingerprints captured before removing the embedded records:
every decoded JSON value, record order, repeated record and metadata field is
preserved.

In a JSON fixture, a SAM string is represented as exactly
`{"sid_member": "<member>", "sam_sha256": "<digest>"}`. `read_fixture` reads
plain or compressed JSON and restores these references; `restore_records`
does the same for an already parsed value. Both leave ordinary JSON values
unchanged. A missing member or a digest absent from that member fails rather
than searching other members. The digest is SHA-256 of the ASCII SAM text,
without a trailing newline. It is distinct from Osteosarc's lossless BAM
record identity: `records` first verifies the exported BAM record identities
and multiplicities against the packaged manifest, then builds the SAM lookup.

## Making the reads again

osteosarc makes Isovar's bundle from a recipe: Isovar's `isovar/` members and
the published T2 small-variant members of a shared bundle's recipe, with the
sources they use and the complete small-variant target inventory
(`sid_data.recipe`). Every member pins its exact records by checksum.

```sh
# From openvax-v1's records (downloaded once, 28 MB; then offline).
python -m isovar.sid_data build --output /tmp/sid-reads
# Every record acquired again from Sid's original BAMs: slow, needs network access
# and the osteosarc metadata snapshot the recipe names (2026-09-25).
python -m isovar.sid_data build --output /tmp/sid-reads --sid
# A new build compared with the packaged bundle; exits 1 if they differ. CI runs this.
python -m isovar.sid_data check
```

Only a cache that holds that snapshot can use `--sid`: a new snapshot has a new
ID, and the recipe pins each BAM to the old one
([iskandr/osteosarc#89](https://github.com/iskandr/osteosarc/issues/89)).

`--from` takes Isovar's members from another bundle: a newer published one, or
a folder such as `isovar/data/sid-reads` itself. `check` compares the members,
their records and the sources' exported headers, not file bytes (the toolchain
versions in a manifest can differ). A bundle carved from another archives its
sources' exported headers as their originals, and records the shared bundle's
source BAMs rather than the original acquisitions; two of openvax-v1's 39
sources thus lose `GO:none` from their archived header
([iskandr/osteosarc#90](https://github.com/iskandr/osteosarc/issues/90)). `build` prints the new manifest's SHA-256. To package a new
bundle, put the folder in place of `isovar/data/sid-reads` and pin that hash as
`PACKAGED_MANIFEST_SHA256`.

Reads are still selected in osteosarc, not here. Tests may consume additional
existing shared members without reselecting their reads. To change the selection,
add it to osteosarc's openvax recipe and release a new bundle; then
build Isovar's bundle `--from` it. Check a local fixture against a bundle with:

```sh
osteosarc test-data check isovar/data/sid-reads fixtures.json
```

`fixtures.json` maps member names to files, or to `{"json": path, "pointer": "/records"}`
for reads embedded in JSON. Decode reference-based fixtures with `read_fixture`
before passing their JSON to this external checker, as the Isovar test does.

The regional audit builders under `tests/data` still acquire larger research
inputs through osteosarc (`sid_data.open_dataset`, `extract_regions` and
`fetch_metadata`), configured with `ISOVAR_SID_SNAPSHOT` and, optionally,
`ISOVAR_SID_CACHE`. Their outputs stay outside the package.
