# GDC viral sequence recovery validation

The installed 1.56.0 wheel was exercised through `python -m isovar viral` on the
actual GDC GRCh38.d1.vd1 FASTA. The archive was acquired from GDC and verified
against published MD5 `3ffbcfe2d05d43206f57f81ebb251dc9`; its extracted FASTA
has SHA256 `fbde121e898f00cc064eaafe5acc58940d7fe069ad14d34acd40c5d675e6b69b`.
The default command checked all 2,779 BAM/FASTA dictionary entries and the DNA
MD5 of every one of the 200 explicitly listed viral contigs. No GTF, gene list
or mutation was supplied.

`examples/benchmark_viral_sequences.py` creates synthetic RNA reads on this
reference, with independently specified overlapping paths, a sequence allele,
a splice, opposite stored orientations, human background and an ambiguous
host/viral placement. This is a functional check, not a biological sample or
evidence of a particular virus in a tumor.

| Evidence | Expected | Recovered |
|---|---:|---:|
| HPV16 paths (reference sequence and allele) | 2 × 400 nt | Exact |
| EBV spliced path | 300 nt | Exact |
| Viral-associated fragments | 87 | 87 |
| Unassigned viral-associated fragments | 0 | 0 |
| Total input alignment records | 108 | 108 |

HPV16 has 74 aligned fragments, including one with an alternative host
placement; 73 have unambiguous viral-contig attribution. EBV has 13 aligned
fragments. Counts conserve original fragments across sequence ambiguity groups,
and the host/viral alternative remains explicitly ambiguous. The FASTA contains
only the three recovered paths. Transcript ends, RNA strand and coding models
remain unassessed. The other panel entries have no eligible alignments in this
synthetic input; the report does not call them absent from biological samples.

The command completed in about 4.5 seconds with the reference index available;
this is a small functional experiment, not a general performance benchmark.
Discovery retained at most eight distinct observations for one contig, with
temporary mirrored views for bidirectional extension. The disk index visits
background once and processes discovery observations one contig at a time.
Retained paths and full-input target witnesses still determine report memory.

The experiment summary in
[`data/gdc-viral-sequence-experiment.json`](data/gdc-viral-sequence-experiment.json)
records acquisition, FASTA/BAM, source and tested-wheel checksums, normalized
commands, exact reconstructed-sequence hashes and per-contig counts. Reproduce
with the supplied local FASTA:

```sh
PYTHONPATH=. python examples/benchmark_viral_sequences.py \
    --genome-fasta GRCh38.d1.vd1.fa --output-dir /new/experiment-directory
```

An offline fixture independently acquired NCBI `NC_007605.1`, with source and
sequence checksums. Its 171,823 bases match GDC's `chrEBV` MD5
`6743bd63b3ff2b5b8985d8933c53290a`. Regression tests recover an exact overlapping
RNA region from it without network access. Other offline cases cover alleles,
splices, late support after discovery limits, per-contig discovery, related
placements, barcode/read-group/mate identity, deduplicated original payloads,
unassemblable reads, filtered records, wrong-reference checksums and CLI export.

This validation also exposed and fixed two shared defects:

- [#509](https://github.com/openvax/isovar/issues/509): missing declared SA pieces
  or mapped mates could incorrectly retain unique attribution. Their evidence
  is retained, with an explicit unresolved-placement reason.
- [#510](https://github.com/openvax/isovar/issues/510): the target-oriented index
  lacked mirrored views required by 5-prime extension. A competing HPV16 path
  stopped at 250 nt despite reads supporting all 400 nt. Both directions now
  extend while preserving original segment identities and counts.

The release adds no Python dependency and uses PyEnsembl's ordinary genomic-DNA
API without GTF indexing. The repository release environment passes `pip check`.
General viral coding annotation and real multi-sample validation remain in #465.
