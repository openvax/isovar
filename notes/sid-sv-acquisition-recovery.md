# Recovery of the failed Sid 10x batch

All **13 targets / 26 orientations** in the pilot's failed 10x batch now have
acquired inputs and offline reconstruction results. This follow-up preserves the
[original four-platform pilot](sid-sv-four-platform-pilot.md) and its failure
receipts. It is not a complete-catalogue result.

## What failed

The original shared extraction exited unsuccessfully; its receipt did not retain
stderr. That specific failure did not recur. A retry instead reached the
500,000-record acquisition cap after 46 seconds. An independent indexed
`samtools view -c -M` count found **1,044,612 regional alignments** in the same
union of intervals. This is regional context, not junction-supporting evidence.

Splitting the batch while preserving its complete target and interval union
recovered twelve targets below the existing cap. Only `::IGKC` exceeded the cap
individually: its 773,288 records made up 74% of the original batch. Streaming
recovery with a larger recorded budget acquired that target in 23 seconds with
about 1.1 GiB peak RSS. No reads, tags, qualities or requested context were
removed to make acquisition succeed.

## Per-target results

Regional counts below were independently recounted from the pinned recovered
BAMs using each target's full interval union. Both orientations were rerun
offline; their scientific results and detailed reconstruction hashes agree with
the first recovered run.

| Target | Regional alignments | Reconstruction in both orientations | Remaining limits |
|---|---:|---|---|
| SV0096 | 43 | No candidate paths | No supplied reference models |
| TRDV1::TRAC | 23,795 | Splice ambiguous | Extension segments |
| ::CNN2 | 98,966 | No candidate paths | 50,000-record reconstruction cap |
| CNTNAP2::CNTNAP2 | 85 | No candidate paths | None reported |
| EEF1AKMT2::FAM53B | 13,804 | Splice ambiguous | None reported |
| UBE2Q2P6::CSPG4P10 | 858 | Splice ambiguous | None reported |
| SV0271 | 1 | No candidate paths | No supplied reference models |
| CD38-associated SV0406 | 26,403 | No candidate paths | None reported |
| POLR2J2::POLR2J3 | 24,901 | No candidate paths | None reported |
| ZFAND2A::C7orf50 | 25,184 | No candidate paths | None reported |
| CDK11A::CDK11B | 35,301 | Splice ambiguous | Extension segments |
| LRRC37A::LRRC37A2 | 21,983 | Splice ambiguous | None reported |
| ::IGKC | 773,288 | Event-linked candidate status | Records, seed segments, extension segments |

These are reconstruction statuses, not confirmed fusion proteins. Some
splice-ambiguous paths yield no exported ORF. A capped search with no candidate
paths does not establish absent RNA support. Even an event-linked result can
include ORFs whose own junction interpretation is splice ambiguous.

The strongest checked `::IGKC` candidate is the short sequence
`MKTDGAATVRLISTLV`: exact full-interval witnesses occur in **487 read fragments
across 99 cell labels** within the collected input. Its junction is splice
ambiguous, RNA polarity is unresolved, and initiation and translation are
unobserved. These are not counts across all 773,288 input records, nor proof of a
novel SV-derived protein or independent RNA molecules.

The aggregate is **14 no-candidate-path views, 10 splice-ambiguous views and two
event-linked views**, with no acquisition-unassessable views in this batch.

## Fixes and validation

- [osteosarc #110](https://github.com/iskandr/osteosarc/pull/110), released as
  0.15.3, streams ordinary seed context while retaining original records and
  recovery provenance. The 5,669,578-record 10x FTH1/FTL input also acquired
  successfully in 235 seconds at about 4.2 GiB peak RSS. UBC/RPS27A and
  TMSB10/TMSB4X acquired 1,351,544 and 4,204,441 records respectively.
- [osteosarc #111](https://github.com/iskandr/osteosarc/pull/111), released as
  0.15.4, fixes the independently reproduced `Illegal seek` cleanup diagnostic
  when a timeout interrupts a BAM header on a pipe. It preserves the actual
  timeout/error and does not claim to fix slow remote queries.
- Isovar's audit now partitions multi-target subprocess failures/timeouts as
  well as record overflow, preserving the error and exact partition. The CLI
  exposes `--max-acquisition-records` and `--seed-timeout`; budgets stay in the
  receipt. Reattempts use a new audit directory rather than overwriting failures.

The [machine-readable ledger](data/sid-sv-recovery/recovered-batch.json.gz)
accounts for all 26 views, records input/result hashes and limits, and includes
five independent candidate checks. Those checks verify 557 witness intervals,
518 transcript translation comparisons, and 495 original BAM record hashes.
Counts describe comparisons performed, not distinct transcripts or proteins.

Reproduce the ledger with the preserved acquisition directories and local BAMs:

```bash
python -m examples.check_sid_sv_recovered_batch \
  /tmp/sid-sv-audit-compact-pilot /tmp/sid-sv-13-target-retry \
  /tmp/sid-sv-dense-retry /tmp/sid-sv-recovery-confirmation OUTPUT
```

The checker requires every intended target, verifies original BAM/index and
upstream provenance pins, independently counts regional records, compares both
offline reconstruction runs, translates candidates/references with Biopython,
and verifies claimed witnesses against original SAM records and the pinned BAM.
It never treats a missing result as an empty result.

The separate remaining discovery/support limit is tracked in
[Isovar #437](https://github.com/openvax/isovar/issues/437). Acquiring all context
does not by itself make reconstruction or per-cell support counting exhaustive.
Broader ONT retries and the full catalogue remain separate from this completed
13-target batch.
