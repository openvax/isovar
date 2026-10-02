# Sid structural-variant RNA audit

This runner freezes all shared SV nominations and the broader fusion survey,
then reconstructs each retained-flank geometry in both directions for each RNA
processing product. It calls Isovar's existing reconstruction and ORF exporter.
No caller-support, gene-expression or protein-coding filter removes targets.
The audit seeds event-related paths, including compatible junctions and
ordinary-splicing ambiguity. It disables unrelated regional seeds with
`include_regional_candidates=False` (`isovar sv-rna --event-only`), while retaining
the same context reads and overlap assembly. Skipped seeds are still reported.

The pinned October 2, 2026 inventory in `data/` contains:

- 637 shared SV nominations, including 120 unresolved geometries;
- 1,621 fusion calls representing the website's 1,271 surveyed junctions;
- 1,493 distinct retained-flank geometries across both inputs;
- 913 alignment products, each with a scope decision and original metadata;
- 77 candidate tumor RNA products, including unavailable and reprocessed files.

The tumor matrix has **173,866 nomination/product pairs**. These are not
independent mutations or libraries. Geometry grouping shares a query, not an
allele verdict: inserted-sequence and caller differences remain in each original
nomination. RNA-only calls are not promoted to DNA-confirmed SVs.

## Reproduce

Run from the repository root with Isovar's dependencies installed. Inventory
acquisition requires osteosarc 0.15.1; offline tests never fetch metadata or reads.
Create an explicit osteosarc metadata snapshot first. The original run used
`sid-sv-audit-2026-10-02`, snapshot ID
`7526b8ab16a806769d5f1868df9931e8dd45d74c2cec4eef0bfe4e59267fdb27`.

```sh
python -m examples.sid_sv_audit inventory audit --snapshot SNAPSHOT --cache CACHE
python -m examples.sid_sv_audit references audit --cache CACHE
python -m examples.sid_sv_audit headers audit --cache CACHE --workers 4
python -m examples.sid_sv_audit acquire audit --cache CACHE --workers 2
python -m examples.sid_sv_audit run audit --cache CACHE --workers 2
python -m examples.sid_sv_audit report audit --cache CACHE --output report --require-complete
```

Alternatively copy the two pinned inventory files from `data/` into a new audit
directory. Reference preparation uses an already installed Ensembl 115 GTF and
cDNA/ncRNA FASTAs. It saves transcript sequences, exon/CDS intervals, source
checksums and explicit exclusions. Noncanonical CDS models retain their sequence
and exons while their unsupported coding frame is recorded separately.

Acquisition requests indexed breakend windows, transcript exons and linked
original mate/SA records through osteosarc. Original qualities and labels are
preserved. It never falls back to downloading an entire remote BAM. Resource
limits and failed requests are recorded. A successful local query cannot prove
that a partner absent from an incomplete acquisition was absent upstream.
Partner queries default to batches of 64 regions, at most 128 queries per
acquisition batch, with a 30-second timeout each. These are explicit resource
bounds, adjustable with `--partner-batch-size`, `--max-partner-queries` and
`--partner-timeout`. Both timed-out and budget-truncated inputs remain usable
but carry acquisition limitations through the ORF export and report.
Use `acquire --source-id ID --batch-id ID` to distribute or pilot individual
acquisitions without changing the intended matrix. Acquisition stops scheduling
new batches when the cache volume falls below `--min-free-gib` (default 8).
This leaves unattempted work pending; it does not manufacture a terminal outcome.
If a request exceeds its record budget, the runner recursively partitions its
geometries, preserving every requested window and transcript exon. The saved
partition tree accounts for each geometry exactly once and pins each child
receipt. A single geometry still exceeding the budget remains an explicit
failure. Offline reconstruction and reports verify the entire partition tree.

`run` operates on verified local BAMs and pinned references; it needs no network.
Checkpoints bind the input hashes, engine source hashes, version and parameters.
Use a new run directory when changing frozen inputs or parameters. A partial
acquisition or wall-time limit is an explicit outcome, not a negative RNA call.

## Read the results

`coverage.json.gz` includes every intended nomination/product pair. A missing
reconstruction keeps the audit incomplete. Unresolved geometry, incompatible
reference builds and missing indexes are recorded as unassessable. Acquisition
limitations remain separate from the reconstruction status.
`all_pairs_accounted` means every pair has a terminal outcome, including
unassessable or limited outcomes; it does not certify unrestricted read recovery
or biological absence. Reports reject mixed reference snapshots, engine versions
or reconstruction settings and verify detailed result checksums.

`candidate-orfs.json.gz` retains exact nucleotide alternatives, protein sequences,
all contributing paths, full-interval witnesses, cell/UMI and read-lineage
evidence, start-tier ambiguity and limitations. The CSV is a convenient view of
the same event-crossing ORFs. Detailed reconstructions also retain regional-only
hypotheses, annotated-frame translations, searched intervals and original reads.
Support is never summed across products. Native-query orientation is recorded
RNA polarity evidence, not proof of RNA strand or translation. A complete ORF
sequence is a candidate, not demonstrated protein expression or presentation.

Blood RNA is a separate cohort (`--cohort blood_control`). Unknown specimen
identity remains separate (`--cohort unresolved_specimen`). The default tumor
selection includes uncertain RNA path hints for header inspection, without
inventing library/timepoint metadata. Immune-repertoire-only products are
explicitly excluded from genomic SV reconstruction.

Known upstream metadata drift is tracked in
[osteosarc #103](https://github.com/iskandr/osteosarc/issues/103); overlapping
Tempus evidence remains tracked in
[osteosarc #100](https://github.com/iskandr/osteosarc/issues/100).
