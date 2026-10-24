# Dense SV evidence-group budget: original RNA replay (#464)

This benchmark records the effect of `max_evidence_groups` on dense discovery.
It uses original COLO829 short-read alignments for SV3 and HCC1395 ONT direct
RNA alignments for the SUGT1 event. The alignments are **selected-site subsets**
of public runs SRR8615617 and SRR23456709, not whole-transcriptome inputs.
The HCC1395 90,603-alignment subset is produced by
[`align_full_hcc1395_rna.sh`](../reports/hcc1395_colo829/align_full_hcc1395_rna.sh).
Input and model SHA-256 hashes, per-path sequence hashes/support, exact counters,
and raw timings are in [benchmarks.json](data/sv-evidence-groups-464/benchmarks.json).

`--all-aligned-regions` marks a **stress case**: it makes every aligned region
in the selected-site BAM eligible as discovery context. The nominated SUGT1
event and original-read support rule are unchanged. This is not a biological
claim about all HCC1395 SV transcripts.

| Source and scope | Record cap | Group cap | Scanned originals | Eligible records / groups | Selected records / groups | Wall s | Peak RSS MiB | Peak scratch MiB | Paths |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| COLO829 SV3 regions | 50,000 | none | 52,212 | 16,641 / 12,820 | 16,641 / 12,820 | 10.80 | 332 | 21 | 20 |
| COLO829 SV3 regions | 100,000 | none | 52,212 | 16,641 / 12,820 | 16,641 / 12,820 | 9.95 | 331 | 21 | 20 |
| COLO829 SV3 regions | 250,000 | none | 52,212 | 16,641 / 12,820 | 16,641 / 12,820 | 9.45 | 331 | 21 | 20 |
| COLO829 SV3 regions | 50,000 | 5,000 | 52,212 | 16,641 / 12,820 | 6,531 / 5,000 | 5.72 | 223 | 20 | 20 |
| COLO829 SV3 regions | 250,000 | 5,000 | 52,212 | 16,641 / 12,820 | 6,531 / 5,000 | 5.47 | 225 | 20 | 20 |
| COLO829 SV3 regions | 250,000 | 12,000 | 52,212 | 16,641 / 12,820 | 15,563 / 12,000 | 15.58 | 320 | 21 | 20 |
| HCC1395 SUGT1 regions | 50,000 | none | 90,603 | 29 / 23 | 29 / 23 | 16.52 | 201 | 186 | 11 |
| HCC1395 all aligned regions | 50,000 | 5,000 | 90,603 | 81,343 / 81,216 | 5,017 / 5,000 | 57.95 | 2,611 | 219 | 11 |
| HCC1395 all aligned regions | 250,000 | 5,000 | 90,603 | 81,343 / 81,216 | 5,017 / 5,000 | 65.24 | 2,614 | 219 | 11 |

The COLO829 50,000/100,000/250,000 record-only runs produced identical
candidate sequence hashes and complete-input support. The 5,000-group runs
also matched each other across record caps; they retained a different set of
20 bounded paths from the unbounded-group runs. At 12,000 groups, 10 of 20
path IDs matched the unbounded result with identical sequence hashes and
support. The event junction still had three original supporting fragments
at every COLO829 setting. The HCC1395 broad-context 50,000 and 250,000
record-cap runs had identical paths and support, since 5,017 selected records
never reached either record cap. Both counted all 90,603 original alignments.

These are single measurements on a shared host; other tests ran concurrently
with some cases. Timings should be treated as resource observations, not a
reliable speedup estimate. In both regional cases, 50,000 exceeds the eligible
discovery count, so 100,000 and 250,000 cannot test additional records. The
HCC1395 broad case has fewer than 100,000 original alignments, and a
250,000-record cap cannot add input. A broad HCC1395 record-only 50,000 run
was not attempted: the 5,000-group run already peaked near 2.6 GiB and a
larger hypothesis set would make this shared-host comparison unsafe.

The budget controls **which exact allelic/structural groups discovery can
search**, not how many original fragments support a retained path. Different
candidate sequences after a tighter budget are an expected loss of search
scope; an omitted group is never negative RNA evidence. Grouping shares
computation but does not deduplicate molecules or calibrate probabilities.

To repeat a case, run a fresh process for each parameter set so peak RSS is
independent. Supply a local HCC1395 selected-site BAM with the SHA-256 in the
JSON file:

```bash
python -m reports.hcc1395_colo829.benchmark_sv_evidence_groups \
  colo829 --max-records 50000 --max-evidence-groups 5000
python -m reports.hcc1395_colo829.benchmark_sv_evidence_groups \
  hcc1395 --bam /path/to/SRR23456709.all_sites.grch38.bam \
  --max-records 50000 --max-evidence-groups 5000 --all-aligned-regions
```
