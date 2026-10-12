# Complete the original pilot's acquisition-failed SV targets

Follow up [#436](https://github.com/openvax/isovar/issues/436), using the source
date validation from [osteosarc #113](https://github.com/iskandr/osteosarc/pull/113).

1. Derive the exact 22 source/geometry pairs (44 orientations) marked
   unassessable in the frozen 400-view pilot. Retain all original outcomes and
   the historical 13-target recovery as provenance. Reconcile later local
   results and reuse only checksum-verified inputs with the same requested
   regions/source identity; vanished inputs require explicitly fresh acquisition.
2. Use the frozen inventory/transcript-reference snapshots. After upstream
   0.15.5 ships, acquire one source/geometry at a time through its public indexed
   extraction/recovery API. Keep every requested context record, original tag,
   quality and multiplicity. Record seed/partner budgets and acquisition failures;
   timeout/cap/missing source results are unknown, never zero RNA support.
3. Maintain at least 8 GiB free and run each source serially within its cache;
   independent attempts require separate mutable cache metadata. Preserve durable original reads
   and archives; discard only regenerated validation-restore copies after their
   retained archive identities/checksums are verified. Never download whole
   libraries. Retry into separately identified attempts when budgets change.
   Regenerable export directories may also be reclaimed after every member of
   their retained archive is verified against the immutable manifest. Keep
   staged coverage snapshots under content-derived names so pending and
   completed states both remain inspectable (#456). Preserve the producing
   code/run pins when reusing an earlier verified acquisition.
4. Reconstruct both orientations with the released dense support engine and
   frozen settings. Full-input candidate support and bounded discovery are
   distinct. Preserve ambiguous splice paths, incomplete partner recovery and
   unknown DNA origins. Check exported sequence translation and selected
   full-interval witnesses independently from original BAM records.
5. Publish a pair/view completion ledger, RNA/protein/support rows and readable
   DNA-event/transcript anatomy, including unresolved outcomes. Integrate the
   new selected results with the existing neoORF report and save portable pinned
   bundles. The full-catalogue denominator remains unchanged.
6. Version bump, lint, full tests, green CI, merge and clean-master deployment
   are required. Close #436 only when all its original acquisition-failed pairs
   are accounted for with truthful revised outcomes.

The first ONT ::IGKC seed request exceeded its 1,200-second budget. Preserve
that immutable failure under #436. Check the native regional count and retry
in a separate attempt with a declared longer seed budget and an isolated cache,
allowing read-only remote retrieval alongside the earlier attempt. Keep the inventory, intervals, record policy,
scientific engine and reconstruction budget fixed. Imports must verify the
actual producing request and explicitly declared per-pair exception; they must
not rewrite the earlier request to resemble the retry.

The native SAMtools/HTSlib 1.21 IGKC retry also stalled in remote HTTP reading.
Preserve its interrupted failure and read-only stack profiles. A separately
identified attempt uses the already installed pysam SAMtools/HTSlib 1.24 through
a pinned subprocess adapter, retaining the same 3,600-second seed budget and
scientific request. Its actual reader version, adapter/library checksums and
reused header are recorded; this does not attribute the stall to a proven
HTSlib bug.

Before IGKC reconstruction, declare a separate 7,200-second per-view attempt:
its complete regional seed has 2,715,423 records, while the 1,758,048-record
TMSB forward scan used nearly the original 3,600-second budget. Preserve the
completed scale-control result and the original budget. Import checks accept
this explicit IGKC-only wall-time exception with unchanged scientific settings;
all producing requests remain intact. Start another dense scan only when its
estimated scratch space leaves the 8 GiB reserve.

Reject mixed-checkout public CLI invocations before preparing or fetching inputs
(#462). One acquisition-only attempt retained an actual Isovar 1.45.2 producer
label; its complete package profile differs only in the version assignment,
and its acquisition code matches published osteosarc 0.15.5. Preserve that
producer and all original requests, require exact source/region/policy and
BAM/index verification for reuse, and reconstruct with the pinned 1.45.3 engine.

For local export, use verified independent APFS file copies to avoid allocating
another complete uncompressed input copy. Keep the public API's original
source/destination checksum verification and the portable archive format. The
public copying/capacity API follow-up is #461. Validate every archived payload
byte, relocate the complete bundle, and test installed CLI archive restoration
and offline reconstruction with a complete selected view, retaining the full
44-view manifest and all original data.
Before writing the recovery archive, require capacity for its entire payload,
conservative tar/gzip overhead and the 8 GiB reserve.

Completed byte-identical cached receipt copies may use independent APFS clones.
Large completed JSON receipts may also use transparent filesystem compression
after verifying their logical size, every byte and original modification time.
Keep files unchanged when the system does not apply compression. Preserve the
storage operation's checksums and helper code; these operations do not alter
read files, scientific inputs or recorded producing requests.

## Re-plan after the FTH forward timeout

The ONT FTH forward view exceeded its recorded 3,600-second limit. Preserve
that result and stop the unfinished reverse scan in the same attempt. Retain
all 22 acquired subsets and all 40 completed views. Use a fresh publication
directory so the original timeout is never overwritten by a later result.

Declare separate FTH and IGKC reconstruction attempts with 14,400 seconds per
view before starting them. Only the wall-time budget changes: the inventory,
references, reads, scientific engine, bounded discovery and full-input support
counting remain identical. Run the two pairs serially and require enough free
space for the observed scratch demand plus the 8 GiB reserve. Preserve all
earlier budgets, controls and unfinished outcomes in the new publication's
provenance. Generate the final report only from completed selected results.

## Re-plan after the working-space stop

The serial coordinator stopped the long FTH reverse scan when free space fell
to 6.69 GiB. Preserve its stop receipt, log and completed forward result; the
interrupted reverse view has unknown support. Do not overwrite the original
coordinator or logs. Import the completed forward view with the same strict
scientific-input checks, then make a separately labeled interim report from
the 97 completed SV views. List the three missing views explicitly and keep
the final 44-view verification/publication gates intact.

Resume the unfinished reverse and IGKC views only with sufficient observed
working storage, using a separately pinned coordinator and the existing
scientific engine and 14,400-second budgets. Reuse the completed forward
checkpoint rather than repeating it. Final report publication, CI, merge and
release remain pending. Investigations of evidence grouping (#464) are
separate from the frozen recovery's sequences and support counts.

## Publication scope agreed on October 11

Publish the 42 completed views and explicitly retain the two unresolved ONT
IGKC views, as requested by the user. The FTH reverse view completed after the
storage stop. The IGKC forward view exceeded its 14,400-second budget; its
reverse scan was interrupted by the coordinator after that failure. Neither
outcome establishes zero RNA support. Preserve the timeout bytes and derive a
separately identified interruption record from the pinned coordinator stop
receipt. No further dense reconstruction is part of this release.

Use a new publication directory, leaving earlier attempts unchanged. Verification
must recount all 22 acquired inputs, validate all 44 recorded outcomes, and
independently check sequences/support for the 42 completed views. Its ledger
must distinguish accounted, reconstructed, and unresolved views. The expanded
55-event report retains all 100 selected SV outcomes; the original pilot keeps
its 400-view denominator. Export 44 outcomes and 42 executable replay plans,
verify the archive and installed-wheel replay, and track IGKC recovery separately.
