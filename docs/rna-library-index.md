# Reusing original RNA libraries

`prepare_rna_library` creates a reusable index of one original BAM. Dense SV
reconstruction, targeted antigen isoforms and viral RNA reconstruction use this
same template acquisition layer. It preserves every original mate and alternative
placement under its `(read group, query name)` identity. A nearby alignment selects
a complete template; it does not turn the unsequenced insert into observed RNA.

## Python API

```python
import pysam
from isovar import ReadCollector, prepare_rna_library, reconstruct_sv_rna

collector = ReadCollector(min_mapping_quality=1, use_soft_clipped_bases=True)
prepared = prepare_rna_library(
    "tumor.bam", "rna-indexes", reference_name="GRCh38",
    read_collector=collector)

with pysam.AlignmentFile("tumor.bam") as bam, prepared.open() as library:
    results = [
        reconstruct_sv_rna(
            bam, **event_inputs, read_collector=collector,
            dense_support=True, library_index=library)
        for event_inputs in nominated_event_inputs
    ]
```

`nominated_event_inputs` contains the ordinary oriented events, reference models,
source identities and budgets documented in [SV RNA reconstruction](sv-rna.md).
The index does not infer those inputs or change their meaning.

`PreparedRnaLibrary` is an immutable descriptor with detached `identity` and
`statistics` snapshots. Its `open()` method returns a context-managed
`RnaLibraryIndex`; open separate contexts for concurrent readers. Reusing a closed
context, changing the original BAM/index bytes, or changing collector eligibility
is an error. Separate source/library labels remain separate analyses even when
their file bytes happen to be identical.

For custom `read_filter` functions or overridden `alignment_filter_reason`
methods, supply a reproducible `read_filter_identity` when preparing and using the
index. That identity must include relevant configuration and closed-over data.
Assembly limits, soft-clip use and event-specific budgets remain separate from
original-record eligibility.

## CLI and cohort use

```bash
isovar sv-audit --manifest inventory.json --output-dir sv-results \
  --library-cache-dir rna-indexes
```

Dense audit views automatically share an index per source/reference context.
`--library-cache-dir` retains verified indexes for later invocations. Without it,
the audit uses temporary indexes shared across that invocation. Verified resume
does not prepare indexes for already completed views. The cache location is an
operational choice and does not change the result fingerprint or scientific
evidence. Cohort SV tasks use the same audit runner.

The current change reuses reconstruction acquisition and dense support queries.
The separate path-assessment stage still needs inventory-wide sequence routing
to remove its remaining per-event whole-library scans; this remains tracked in
[#593](https://github.com/openvax/isovar/issues/593).

## Scope, identity and resource use

- The original BAM remains authoritative. The persistent index stores virtual
  offsets and conservative alignment bounds, without a second full BAM copy.
- Build identity includes original bytes, collector eligibility, declared
  assembly and implementation/dependency identity. Incomplete builds are never
  published; cached metadata and database hashes are checked before reuse.
- Interval queries can overselect skipped introns or rounded spatial bounds.
  Existing exact CIGAR, sequence, junction and placement checks determine support.
  A soft clip does not expand genomic coverage.
- Excluded records and exact duplicate payloads remain in original-input
  accounting. Identical SAM records do not create additional support. Different
  sequences, qualities, labels or placements remain distinguishable.
- Template iteration has stable `(read group, query name)` order. Bounded
  variant-free discovery now declares that order, instead of depending on a
  temporary collation's physical order. Complete support still considers every
  eligible selected template; discovery limits remain explicit.

Index preparation takes one source-wide indexing pass and a duplicate-payload
pass over repeated template identities. Disk usage scales with original record
count. Subsequent event views query regional candidates and retrieve their full
templates. Ordinary calls without a supplied index use the same implementation
with a private lifetime; SAM/CRAM inputs are converted to a temporary seekable BAM.

The SAM coordinate/record contract is defined by the
[SAM specification](https://samtools.github.io/hts-specs/SAMv1.pdf). Conservative
spatial bounds follow [SQLite R-tree rounding](https://www.sqlite.org/rtree.html).
