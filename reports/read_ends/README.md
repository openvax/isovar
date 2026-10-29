# Original-read annotation reuse

The 1.71.0 receipt `cache_benchmark.json.gz` compares repeated annotation of four unchanged Sid
records from the packaged, source-linked fixture: the adapter-containing Illumina
R2, two ONT tail candidates and the processed PacBio missing-QUAL record.

Run `python examples/benchmark_read_end_inference.py --iterations 200 --output
benchmark.json.gz` from the repository. The receipt records exact exported BAM
hashes, query identities, the profile configuration/hash and driver hash.

For 800 requests, this run measured 0.02439 seconds without reuse and 0.00618
seconds with reuse: four initial misses and 796 hits, retaining 4,377 sequence
bases. These are repeated annotation timings on four selected records, not
whole-BAM throughput or an estimate of adapter prevalence. All derived views
matched, and every original SAM record remained unchanged.

The documented R2 sequence is an explicit motif hypothesis. It does not identify
Sid's actual library preparation. No resulting false fusion is claimed. The
processed PacBio record is a preservation control, not an adapter-positive
sensitivity validation.

## Placement-aware 1.72.0 replay

`internal_policy_cache_benchmark.json.gz` repeats the same original records with
fresh CIGAR/MD origin assessment after cached sequence matching. It measured
0.03404 seconds uncached and 0.01544 seconds cached for 800 requests (4 misses,
796 hits, 4,377 retained bases). All derived views agree and the original records
remain unchanged. These timings are not an end-to-end benchmark and should not
be compared across runs as a calibrated speedup. The receipt pins the new origin
module as well as the profile, matcher, collector and driver.

Internal linker positives and genomic lookalikes are explicitly simulated in
`tests/test_internal_technical_sequences.py`. They exercise noisy matching,
both orientations, supplied order conflicts, missing/malformed MD, competing
reference claims, supplementary records, CIGAR shortcuts, independent-template
counts and routed/full-library parity. The original SV3 witness replay also
retains its existing IDs and evidence counts when no exclusion policy is enabled.
