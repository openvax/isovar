# Sid test reads

Isovar's tests use original reads from the public Sid osteosarcoma data (CC0).
Those reads now come from **openvax-v1**, the test-data bundle osteosarc
publishes for all four OpenVax libraries (Isovar, Varcode, Vaxrank and
Topiary). Isovar no longer ships or builds its own copy.

## What openvax-v1 holds for Isovar

- **345 members, named `isovar/<path>`** after the files Isovar's tests read
  (`<path>` is relative to `tests/data`). They hold exactly the same records.
- **124 whole read files**, 5.8 MB of SAM and BAM:
  - the six-locus RNA;
  - the 49-case vaccine corpus and its primary-only copies;
  - the stress corpus;
  - the PIP5K1A and PacBio comparisons;
  - the long-read fusions;
  - the chimeric examples.
- **221 selections inside files that stay in the repository**: `file#/pointer`
  names reads embedded in a JSON fixture, and `file#READ` names one read of a
  SAM file.
- **Selection and provenance.** osteosarc owns read selection and records the
  source, snapshot and reason for every template. The bundle also covers all
  187 osteosarc site variants and the other libraries' fixtures. See
  [osteosarc's shared test data](https://iskandr.github.io/osteosarc/test-data/#shared-test-data-openvax-v1)
  and [iskandr/osteosarc#56](https://github.com/iskandr/osteosarc/issues/56).

These are test selections, not a random sample or a whole locus, and they
can't estimate sample coverage or VAF. Historical allele definitions and
protein expectations are unchanged.

## Using the reads

The first time a test needs the reads, `isovar.sid_data` downloads openvax-v1
(28 MB), checks it against the pinned manifest hash, and exports Isovar's read
files in their original formats. The cache then holds about 117 MB. Later runs
are offline:

```python
from isovar import sid_data

path = sid_data.path("osteosarc/bulk_star_t0.sam.gz")   # a local file
sid_data.export("fusions/corpus/TPST1--CRCP-T1.input.json.gz#/original_records", "reads.sam")
sid_data.members()                                      # every Isovar member
```

```sh
python -m isovar.sid_data list
python -m isovar.sid_data path osteosarc/expansion/corpus/28-NTF3-chr12-5494381-53f498a544883d51.bam
python -m isovar.sid_data export chimeric/osteosarc-ont.sam --output /tmp/chimeric.sam
```

The files live in osteosarc's cache: `OSTEOSARC_CACHE`, or the shared OpenVax
cache (`OPENVAX_DATA_CACHE`, or the platform's `openvax` cache directory). CI
caches that folder.
- The exports sit in the cache's `isovar/` folder. They are read-only, and
  checked against their recorded digests the first time each process uses them.
  A damaged export raises an error that names the folder to delete.
- Without network on the first run, the first test to need the reads fails
  with a hint, and the rest fail at once.

An exported file holds exactly the member's original records, but it is
coordinate-sorted and carries the source's full header. So its bytes, and the
order of records at the same position, differ from the old copies. Tests
compare records, not file checksums. `tests/test_sid_data.py` checks every
embedded or selected read against the bundle, so the JSON fixtures can't drift
from it.

## Changing the reads

Reads are selected in osteosarc, not here. To add or change a test's reads, add
them to osteosarc's openvax recipe and release a new bundle; Isovar then moves
to that bundle name. Check a local fixture against a bundle with:

```sh
osteosarc test-data check openvax-v1 fixtures.json
```

`fixtures.json` maps member names to files, or to `{"json": path, "pointer": "/records"}`
for reads embedded in JSON.

The regional audit builders under `tests/data` still acquire larger research
inputs through osteosarc (`sid_data.open_dataset`, `extract_regions` and
`fetch_metadata`), configured with `ISOVAR_SID_SNAPSHOT` and, optionally,
`ISOVAR_SID_CACHE`. Their outputs stay outside the package.
