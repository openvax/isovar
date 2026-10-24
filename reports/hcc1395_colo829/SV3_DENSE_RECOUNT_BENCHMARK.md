# COLO829 SV3 dense recount benchmark

Isovar #562 changes how complete original-read templates are traversed during
dense SV reconstruction. It first makes a private query-name-collated BAM,
then reads each eligible, deduplicated template sequentially. The original BAM
and SAM records are unchanged. The [SAM specification](https://samtools.github.io/hts-specs/SAMv1.pdf)
defines paired, secondary, and supplementary records; [samtools collate](https://www.htslib.org/doc/samtools-collate.html)
guarantees that records sharing a query name are contiguous. Isovar still
separates read groups within that name and counts distinct original fragments.

## Exact SV3 replay

Both implementations replayed the pinned 52,212-record
`SRR8615617.cnih3_source_candidates.grch38.bam` with the same Ensembl 87 CNIH3
models and `max_records=20,000`. Each retained path
`653d5d9c015d3e590fb149b3` with three original split-junction fragments and
three local exact fragments, 50,489 eligible segments, 16,641 selected records,
and `mate_location_unavailable` plus `path_limit` limitations.

| Implementation | Complete dense replay wall time |
|---|---:|
| Isovar 1.59.18, random paired-template seeks | 15.18 s |
| Name-collated sequential template stream | 9.23 s |

The committed original-read regression checks the exact join counts. Timing is
one local run per implementation, not a cross-machine performance guarantee.

## Source-scale paired-load proxy

For a repeatable load test, we wrote 32 copies of the pinned BAM to a private
scratch BAM, prefixing every QNAME with a different `copyNN-` label. The result
has 1,670,784 records and is 155,131,289 bytes. The copies exercise pair and
alternative-placement grouping at approximately the 1.7-million-record target
screen scale; **they are not additional biological molecules or SV3 evidence**.

| Implementation | Spool/index | One complete template traversal | Total |
|---|---:|---:|---:|
| Isovar 1.59.18 | 37.22 s | 96.97 s | 134.19 s |
| Name-collated stream | 53.75 s | 2.62 s | 56.37 s |

Both traversals yielded 863,232 templates, 1,665,024 eligible records,
1,615,648 eligible segments, and 5,760 excluded unmapped records. The larger
one-time collation cost is recovered during the first complete traversal;
dense reconstruction traverses templates again for support. On this load
proxy, the new complete dense SV3 replay finished in 90.86 seconds and found
96 direct and 96 local exact fragments, exactly 32 times the three original
witnesses. It retained a 20-path search limit and a discovery record limit.
These replicated counts must never appear in a biological SV3 report.

Run the pinned benchmark with
`python -m reports.hcc1395_colo829.benchmark_sv3_dense_recount`.
Pass `--templates-only` to time just indexing and traversal, or `--bam PATH`
to use a larger private input. The load proxy can be regenerated with:

```python
from pathlib import Path
import pysam

source = Path("tests/data/colo829_whole_run/SRR8615617.cnih3_source_candidates.grch38.bam")
target = Path("/private/tmp/isovar-562-sv3-paired-32x.bam")
with pysam.AlignmentFile(source, "rb") as original:
    header = original.header.to_dict()
header.setdefault("HD", {})["SO"] = "unsorted"
with pysam.AlignmentFile(target, "wb", header=header, format_options=[b"level=1"]) as output:
    for copy in range(32):
        with pysam.AlignmentFile(source, "rb") as original:
            for read in original.fetch(until_eof=True):
                read.query_name = f"copy{copy:02d}-{read.query_name}"
                output.write(read)
```

This performance change does not resolve the mature CNIH3 exon path or
reading frame. The original five exact-window spots and three placed joins
remain the biological evidence; [#563](https://github.com/openvax/isovar/issues/563)
tracks the missing connection.
