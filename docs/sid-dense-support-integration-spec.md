# Dense support integration for the Sid audit

Apply the released dense discovery/support engine to the existing offline Sid
audit without replacing its original bounded pilot outcomes.

1. Merge the released engine into the audit branch, enable dense support in its
   pinned reconstruction parameters, and retain discovery/support scope in
   summaries, exported candidates and coverage reports. Scratch storage is local
   and disposable; checkpoints must still bind inputs, code and settings.
2. Add regressions through the acquisition/reconstruction/report adapter for
   late event discovery, complete support beyond the discovery budget, explicit
   acquisition incompleteness, and request drift.
3. Reuse verified original pilot acquisitions in a fresh result directory and
   rerun affected dense inputs in both orientations. Independently verify
   original witnesses and per-cell counts, retain input/engine/result pins, and
   publish measured results alongside the original pilot.
4. Run lint and the full suite, update PR #433 and verify its CI. The user has
   narrowed shipment to the focused subsets in `sid-priority-subsets-spec.md`;
   the full catalogue remains pending with its denominator unchanged.

Full-input support means all supplied acquired records were scanned for retained
candidate hypotheses. It does not remove acquisition limits, guarantee complete
hypothesis discovery, or establish initiation, translation or biological absence.
