# Original-read annotation reuse

`cache_benchmark.json.gz` compares repeated annotation of four unchanged Sid
records from the packaged, source-linked fixture: the adapter-containing Illumina
R2, two ONT tail candidates and the processed PacBio missing-QUAL record.

Run `python examples/benchmark_read_end_inference.py --iterations 200 --output
benchmark.json.gz` from the repository. The receipt records exact exported BAM
hashes, query identities, the profile configuration/hash and driver hash.

For 800 requests, this run measured 0.02274 seconds without reuse and 0.00608
seconds with reuse: four initial misses and 796 hits, retaining 4,377 sequence
bases. These are repeated annotation timings on four selected records, not
whole-BAM throughput or an estimate of adapter prevalence. All derived views
matched, and every original SAM record remained unchanged.

The documented R2 sequence is an explicit motif hypothesis. It does not identify
Sid's actual library preparation. No resulting false fusion is claimed. The
processed PacBio record is a preservation control, not an adapter-positive
sensitivity validation.
