# Preserve splice-ambiguous event hypotheses in dense discovery

Issue #441 is reproduced in the Sid ::IGKC replay: `_join_class` deliberately
defers splice-ambiguous event joins as `lazy`, but dense selection ranks them
below unrelated breakpoint clips. Promote both placed event classes before clip
context; keep existing biological classification and ambiguity unchanged.

Regress with abundant breakpoint clips and a repeated event-compatible ordinary
splice behind them. A two-record discovery budget must retain the splice
hypothesis, and support must count every eligible original fragment and cell.
Run lint, the full suite and GitHub CI, then merge/deploy the patch before the
corrected audit replay. Keep the initial Sid replay as failure evidence.
