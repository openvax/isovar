# High-depth RNA performance (#227)

PR: [#230](https://github.com/openvax/isovar/pull/230). No production defaults,
read filtering, genetic-code rules, or public APIs change. Version: 1.8.2.

## Full-depth benchmarks

These are original checksum-verified regional alignments, not the bounded test
fixtures. Baseline production code is the unchanged `v1.8.1` tag. Each row is a
fresh-process run on the same shared macOS/ARM64 workstation, Python 3.12.6.
Wall times are observations under variable load, not portable latency promises.
CPU time is also recorded for final runs. Reference setup, checksum verification,
diagnostic fingerprinting and report serialization are outside production timing.
Peak RSS is the process high-water mark, not additional memory used by a stage.

| Original source / locus | 1.8.1 production | New production | New independent checks | Ranked windows |
| --- | ---: | ---: | ---: | ---: |
| ONT T3 deduplicated / MT-ND5 | 1,386.0 s | 69.2 s | 8.9 s | 2,854 |
| ONT T1 tagged, primary-only / MT-ND5 | 165.3 s | 14.6 s | 2.6 s | 2,006 |
| ONT T1 tagged / DYNC1H1 chr14:101980529 | 6.59 s | 1.50 s | 0.019 s | 241 |
| ONT T2 tagged / H1-2 deletion | 39.87 s | 3.42 s | 0.051 s | 384 |
| Bulk 2024 / DYNC1H1 chr14:102030200 | 0.064 s | 0.092 s | 0.003 s | 112 |
| ONT T3 tagged / MT-ND5 | Not measured to completion | 113.1 s | 24.5 s | 4,351 |

All five completed baseline comparisons preserve the exact ordered RNA
fingerprint (candidate sequences, all original read fields and coverage) and
the complete ranked protein/check/support-name fingerprint. The small bulk
control remains subsecond; its timing does not establish a speedup.
The T3 tagged result is independently validated but is not advertised as a
completed whole-pipeline baseline comparison.

The original 180-second profiles attribute 162–165 seconds to exhaustive
support matching. The independent checker is timed/profiled separately;
the previous T1-primary validation timeout was caused by sharing the audit
budget with slow reconstruction, not a demonstrated translation error.
Intermediate 120-second failures are retained in the benchmark evidence.

Memory remains substantial: T3 deduplicated peaks near 5.2 GB and T3 tagged
near 9.8 GB. Those peaks have already been reached during read collection,
before the support index exists. [#229](https://github.com/openvax/isovar/issues/229)
tracks that separate allocation problem. This PR does not claim to fix it.

## Audit evidence and reproduction

The checked-in [results manifest](results/manifest.json) pins
[benchmark measurements/profiles](results/benchmarks.json) and
[lossless rerun results](results/reruns.json.gz). All **10** interrupted
mitochondrial modes complete: **28,890** pre-cap ranked windows pass independent
protein validation. All **215 nuclear rows / 430 modes** retain their complete
previous results, and counts are unchanged for all 220 rerun rows.
The default top windows are 26–28 amino acids; every mode also has 49-aa
alternatives. That is the unchanged balanced ranking, not full-length protein
reconstruction. There are 64,195 checked translation/transcript contributions;
these overlapping windows are not independent biological observations.

The follow-up reruns all 44 loci for each of the five previously interrupted
RNA products in both default-input and primary-only modes. Its explicit audit
budget is **600 seconds per locus/mode**, versus the historical report's 120.
This larger checking budget is not itself a performance improvement; the
measured production times above are separate evidence.

The packager refuses changed counts, changed completed nuclear results,
unfinished modes, or failed independent protein validation. It retains every
ranked mitochondrial window and exact supporting template-name set, not hashes
alone. Dense membership bitsets avoid hundreds of MB of repeated read names:
`performance_report.supporting_read_names(row, protein)` decodes and verifies
the lossless set against its SHA256. The parent `audit/` matrix stays a historical
1.8.1 artifact; this follow-up does not retroactively rewrite its outcomes.

With the acquisition cache described in the parent README:

```sh
python -m tests.data.osteosarc.expansion.benchmark \
  --destination /path/to/verified-cache \
  --source-id 62a9d3e67e1b564e \
  --output /path/to/new-benchmark-directory --timeout 600
```

Add `--profile` for separate production and validation profiles, or
`--mode primary_only` for the verified flag-filtered derivative. Use
`--variant-id` for a nuclear control. To benchmark 1.8.1, run the same script
with `PYTHONPATH` and the working directory pointing to a clean detached
`v1.8.1` worktree; its identity records the imported production code.

Run `tests.data.osteosarc.expansion.runner` with the same verified cache,
`reference-GRCh38`, `--inventory inventory-GRCh38-validated.json`, and
`--time-limit 600`, once for each source ID:

```
53f498a544883d51  ONT T1 tagged
1120a096937e29e1  ONT T2 tagged
2bc1fc291308debb  ONT T2 deduplicated
62a9d3e67e1b564e  ONT T3 deduplicated
58c4d68e5f7e55b5  ONT T3 tagged
```

Use a fresh output directory. `performance_report --help` describes the
packaging/comparison CLI. Checkpoint hashes, source identities and benchmark
artifact checksums are retained. Production source hashes correspond to commit
`615b453`; baseline production hashes correspond to `v1.8.1`.
Do not run with `python -O`: independent
validation assertions must remain enabled.

These are technical RNA products, not independent patients or molecule/cell
counts. Mitochondrial translation remains conditional on origin: NUMTs, DNA
heteroplasmy and malignant-cell fraction are not resolved by this optimization.

## Spec

Profile the checksum-verified full regional osteosarc sources, not the bounded
offline fixtures. Measure production reconstruction and independent protein
validation separately, including wall time, peak process RSS, and Python
function profiles. Keep timeouts explicit and retain completed evidence.

Optimize measured hotspots without changing the candidate set, original read
objects/multiplicity, coverage, translations, or protein ordering. Compatibility
must remain anchored at the allele and non-transitive: a short shared core does
not connect conflicting longer haplotypes. No downsampling, changed thresholds,
reference padding, or new biological assumptions.

Verify against the exhaustive matcher using randomized and adversarial cases,
and compare complete outputs on original reads. Include empty flanks, deletions,
duplicate names, conflicting branches, coverage trimming, and adaptive lengths.
Benchmark full-depth mitochondrial sources and smaller nuclear controls. Rerun
both modes of the five interrupted #218 mitochondrial cells, with independently
checked proteins and provenance-bearing results. Document remaining limits.

Run lint, the complete offline test suite, and GitHub CI; include a patch version
bump and link the PR to #227 and #218.

## Profile-driven refinement

The first indexed run moved the bottleneck into per-read coverage updates.
Accumulate exact coverage from the disjoint original-read groups while matching,
and carry an independent copy of coverage through interval trimming. Preserve
the per-read interval calculation for standalone `VariantSequence` objects.
Disabled logging must not format coverage arrays eagerly.

The remaining collection-stage memory peak is tracked separately in #229.
Some exploratory runs overlapped other workstation jobs; retain timed-out
measurements and report both process CPU time (where recorded) and wall time.
