# Sid four-platform SV RNA pilot

This is the completed pilot, **not the full-catalogue result**. The pilot covers
50 retained-flank geometries (71 original nominations), four RNA alignment
products and both orientations: **400 intended views**. Every view has an
outcome. A second offline run with the reviewed checkpoint and compact-provenance
implementation reproduced all 400 scientific results exactly; the full candidate
JSON and CSV are byte-identical. The published pins describe this second run.
There are 356 reconstructed views and 44 unassessable views after
bounded acquisition failures. No failed acquisition counts as negative RNA
evidence.

| RNA product | No candidate paths | Splice-ambiguous candidates | Event-linked candidates | Unassessable | Candidate sequences |
|---|---:|---:|---:|---:|---:|
| PacBio T1 (`8275c5f21164`) | 74 | 14 | 12 | 0 | 144 |
| ONT T2 tagged (`58f81602ec2f`) | 42 | 34 | 12 | 12 | 544 |
| Bulk RNA T2 SARC0277 (`40a3ceff825d`) | 42 | 42 | 16 | 0 | 350 |
| 10x T2 Kamil (`6e9147690525`) | 28 | 30 | 10 | 32 | 217 |
| Total | 186 | 120 | 50 | 44 | 1,255 |

The candidate count includes alternative starts, frames and nucleotide sequences
within products. It is not a count of independent SVs, translated proteins or
neoantigens. Support is never summed across products. The pilot takes the first
50 sorted geometry IDs; it is not a random sample or a catalogue sensitivity
estimate.

## Independent checks

[check_sid_sv_pilot.py](../examples/check_sid_sv_pilot.py) independently checks the
highest full-interval-support candidate from each product. It reads original SAM
text, restores sequencing orientation and hard-clip coordinates, verifies every
claimed witness interval, and recounts distinct read-group/query-name fragments,
segments, cell labels and missing qualities. Biopython translates the candidate
and available pinned transcript sequences independently of Isovar's translation
code. The script follows the [SAM specification](https://samtools.github.io/hts-specs/SAMv1.pdf)
and [NCBI standard genetic code](https://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi#SG1).

| Product | Nomination | Exact full-ORF fragments | Read segments | Labelled cells | Segments without qualities |
|---|---|---:|---:|---:|---:|
| Bulk T2 | `fusion-call:1202`, FCGR3B::FCGR3A | 246 | 288 | 0 | 0 |
| ONT T2 | `fusion-call:0932`, ::TRGC2 | 58 | 58 | 55 | 0 |
| 10x T2 | `fusion-call:0414`, ::FAM189A1 | 20 | 20 | 18 | 0 |
| PacBio T1 | `fusion-call:0932`, ::TRGC2 | 16 | 16 | 16 | 16 |

The two TRGC2 representatives have different candidate sequences. The bulk
representative is a short ORF with UTR/noncoding reference interpretations;
the 10x representative is eight amino acids and splice ambiguous. The long-read
representatives have unresolved RNA polarity and sequence-only start evidence.
These examples demonstrate sequence reconstruction and evidence attribution;
they do not establish tumor-specific rearrangements, initiation or translation.
The four checks cover 679 witness intervals and seven transcript comparisons.

The existing pinned TPST1::CRCP regression independently expects
`MPSRRRGGSSLIPGRMPILRFSVTTRYFSY`, a stop codon and 15 original PacBio witnesses
without base qualities. It also verifies that duplicate regional occurrences do
not inflate event support.

## Limits and complete-catalogue denominator

Of the 356 reconstructed views, 248 use explicitly incomplete recovered input.
The 44 unassessable views comprise ONT seed timeouts (five geometries), a 10x
acquisition subprocess failure (13 geometries), and single-geometry record caps
(one ONT and three 10x geometries). Each geometry has two views. Partitioning
oversized requests preserved the other targets and every requested context
region. A remaining pipe-close diagnostic is tracked in
[osteosarc #106](https://github.com/iskandr/osteosarc/issues/106); the subprocess
failure is not assumed to have the same cause.

The full frozen inventory has 2,258 nominations, 1,493 geometries and 77 candidate
tumor RNA processing products: **173,866 nomination/product pairs**. Header
inspection found 44 usable GRCh38 products, five GRCh37 products without matching
audit references, and 28 products without a listed index. The latter groups stay
explicitly unassessable. There are 120 unresolved nominations. Blood controls
remain a separate cohort. Reprocessed products and aliases are not independent
biological observations.

The pilot accounts for 284 nomination/product pairs; 93,788 remain pending
reconstruction, 9,240 have unresolved geometry, and 70,554 are unassessable from
source availability/reference checks. Consequently, the full audit's
`all_pairs_accounted` remains **false**. Full execution and its results report
are required before merging and deploying PR #433.

## Storage measured before the full run

The pilot cache also includes earlier acquisition benchmarks, making it a
conservative starting point rather than an exact per-target cost model.

- Original cache: about 9.5 GiB; compacted cache: 2.6 GiB (2,713,271,878 logical
  bytes after counting hard links once).
- Query-name command inputs: 461 copies totaling 2.26 GB; only 15 distinct lists,
  totaling 68.7 MB.
- All 1,039 original BAM/index hashes remained identical after compacting 546
  receipts. Compaction took 131 seconds and did not reacquire reads.
- Audit acquisition/header metadata: about 2.20 GB before; 7.77 MB with pinned
  upstream receipts and shared assets. Detailed reconstruction/results add
  approximately 158 MB, plus references and the published report.

[osteosarc #108](https://github.com/iskandr/osteosarc/pull/108), released as
**0.15.2**, fixes the duplication while preserving semantic cache identities and
expanded in-memory receipts. The audit now pins full upstream provenance instead
of embedding it again. Missing or changed provenance prevents offline reuse.

Scaling approximately 2.9 GB by `1493/50 × 44/4` gives approximately **0.95 TB**.
This is an extrapolation, not a measured requirement: library depth, context
sizes, failed acquisitions, repeated regions and retries vary. Provision
**2 TB free working storage**, preferably SSD, with a free-space guard. A 500 GB
allocation is too tight. Source BAMs total 851 GiB across the 44 ready products,
but acquisition remains indexed and regional, rather than downloading them all.

## Reproduction and artifacts

Use the [runner instructions](../examples/sid_sv_audit/README.md), the pinned
inventory and Ensembl 115 references. The pilot root batch is
`e11a46eac1e8c5e207fcf974`. Run acquisition once per listed source with
`--batch-id`, then run reconstruction offline. Use a new output directory when
changing engine code or request parameters.

The [machine-readable pilot coverage](data/sid-sv-pilot/coverage.json.gz) contains
all 400 observed views, complete-catalogue counts and input/engine pins. The
[compressed candidate CSV](data/sid-sv-pilot/candidate-orfs.csv.gz) lists all 1,255
candidate rows. The [independent check record](data/sid-sv-pilot/independent-checks.json)
contains the verified counts, candidate identities and original-record hashes.
Full candidate JSON and result pin lists are retained with the local audit;
the published artifact manifest records their checksums.

```sh
python examples/check_sid_sv_pilot.py AUDIT REPORT independent-checks.json
```

Code review found and fixed interrupted-checkpoint publication and alarm timing
([Isovar #435](https://github.com/openvax/isovar/issues/435)). The reviewed changes
pass lint and all **5,707 tests**, with **97% coverage**. Upstream 0.15.2 passes
403 tests and all supported Python/consumer CI jobs.
