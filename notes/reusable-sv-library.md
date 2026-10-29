# Reusable original-RNA library evidence (#593)

Prepare one immutable, checksum/policy/implementation-pinned index per original
RNA library. Index virtual offsets in the original seekable BAM, template keys
(read group, query name, segment), exact duplicate records, eligibility and
conservative genomic spans. Preserve all original records for assessment;
eligibility applies to reconstruction, not to a relabelled subset of the source.
Use a spatial index to find complete templates, then retain the existing exact
CIGAR/clip/event tests. Distant mates and alternative/supplementary placements
must follow their template even when they are outside every queried interval.

Expose preparation/opening as a typed library API and an optional persistent
cache in the existing SV CLI and audit runner. Ephemeral and persistent use must
share one implementation. Publish a validated index atomically; use read-only
connections afterward. Keep per-event discovery tables in a separate temporary
database. Verify source bytes once when opening a library context, and detect
changes during its lifetime. Never reuse an incomplete or incompatible cache.

Audit reconstruction and assessment become explicit phases. First reconstruct
every view using the shared library index. Route original sequences to the union
of retained paths in one whole-library pass, including unmapped and remote
sequence-only matches. Store routed offsets/complete-template retrieval on disk;
do not load the whole library into Python dictionaries. Shared fixed-length seed
matching can use NumPy's exact byte-array search without changing exact-match
semantics. Per-event exon links, alternative placements and statistics remain
specific to that event's reference models. No-path events reuse source counts
and indexed exon queries without loading unrelated sequences.

Keep source libraries/cultures separate, original exclusion and duplicate
accounting, deterministic representative selection, budgets, all witness IDs,
and distinct discovery/junction/exon/frame/complete-sequence states. Source scope
remains the complete original input, with index provenance reported separately.

Verify against the exhaustive implementation on original-RNA fixtures and
adversarial cases: read-group/name collisions, duplicate originals, distant
supplementary pieces, soft clips, unplaced mates, remote sequence-only matches,
changed bytes/policy/code, interrupted preparation and concurrent readers.
Benchmark whole HCC1395 direct RNA and MAS separately; report one preparation,
event-query scaling and peak memory, then apply the shared workflow to the pinned
inventory. Preserve explicit pending states until those analyses finish.

Coordinate/identity definitions follow SAM v1.6:
https://samtools.github.io/hts-specs/SAMv1.pdf . Spatial overlap is only a
conservative acquisition filter; the existing biological interpretation stays
in the reconstruction and assessment modules.
# Shared consumers

The old private spool is also used by antigen and viral RNA reconstruction.
Replace it once with a general `rna_library` index, not an SV-only duplicate.
Keep per-analysis work in private contexts and original eligibility/source identity
in the shared immutable index. Bounded variant-free discovery must declare a
stable template-identity order instead of depending on a temporary collation's
physical order; complete evidence counts remain unchanged. Verify the existing
SV, antigen, viral, SAM/BAM and benchmark consumers together.

## Follow-on: inventory-wide path evidence

The index PR addresses repeated reconstruction preparation, but must not close
#593 while `assess_sv_rna_paths` still scans the BAM for each view. Split the
audit into planned jobs, reconstruction, one source-wide sequence-routing pass,
and per-job assessment/reconciliation. Keep native snapshots/receipts immutable.

Route the union of exact path seeds across all original BAM records and raw ends
once per source. Store matching template/job identities and raw sequence witnesses
on disk. For each view, retrieve complete matching templates plus every template
overlapping that view's reference exons, without collector eligibility filtering.
Use a shared record-to-evidence helper for indexed and direct assessment. This
preserves all alternative placements, exon-link label hashes, raw/unmapped ends,
original record/fragment denominators and source-specific reference exon counts.
Pathless views need indexed exon evidence and original-library counts, not a full
sequence scan. No candidate seed match is itself counted as exact path support.

Validate entire assessment outputs against the current exhaustive collector for
multiple events, empty paths, repeated shared sequences, alternative/supplementary
placements, excluded/no-quality records, raw ends, reverse strands, duplicate
physical records and varying thresholds. Only then benchmark actual inventories
and close #593. A compiled exact seed-routing implementation can replace Python
word iteration if profiling justifies it; preserve literal sequence behavior.

## Inventory runner phases and recovery

Plan and validate every selected job before computation; handle verified completed
views and explicit unavailable outcomes first. Process source/reference groups
serially. Reconstruct into checksum-pinned stage artifacts, then route the union
of retained sequences once and assess each view through the common collector.
Keep staged reconstructions on interruption and verify them on resume; remove
only this job's stage files after its final snapshot and receipt are written.
Preserve original inventory order in the final report. Do not retain every full
reconstruction or raw FASTQ sequence in Python memory across the inventory.

Original-order selected record iteration must include all physical records;
sequence/placement deduplication remains in the assessment layer. Raw-only
sequence matches also select their BAM template's distant alternative placements.
All reference-exon alignments, including ones unrelated to paths, remain available
so exon-label hashes, denominators and complete assessment output stay identical.
# Assessment reuse review

Pin each collected view's reference models and exon threshold, and each routed
view's path sequences and match threshold. Validate those inputs when evidence
is reused. This prevents a direct API caller from silently assessing a new path
or shorter match using an incomplete routed subset. Keep legacy caller-created
evidence dictionaries compatible, and preserve native assessment output.

