# Stream complete SV templates in dense support (#562)

## Goal

Finish complete-input recount on paired COLO829 SV3 candidates without a
random BGZF seek for every original read pair. Keep the existing fragment,
read, UMI, cell, quality, and alternative-placement semantics.

## Method

1. Collate the eligible source BAM by query name once in private scratch.
   Samtools collate groups all records of one QNAME; read-group keys still
   separate distinct templates. Keep every SAM record, including secondary
   and supplementary records, and retain the existing exact-payload duplicate
   check. See the [SAM specification](https://samtools.github.io/hts-specs/SAMv1.pdf)
   and [samtools collate manual](https://www.htslib.org/doc/samtools-collate.html).
2. Scan eligible, deduplicated offsets in collated physical order. Skip
   excluded records sequentially and emit each template once. Keep random
   access only for the small set of cited templates after recount.
3. Compare support and original-record metadata with the current behavior on
   paired, supplementary, duplicate, and cross-read-group fixtures. Benchmark
   a source-sized paired input and replay the pinned SV3 candidate BAM.

## Boundaries

This is an exact recount performance change. It cannot establish SV3's mature
exon path or reading frame without new original RNA linkage evidence (#563).
