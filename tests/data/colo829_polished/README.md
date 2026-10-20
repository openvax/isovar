# COLO829 polished representative allele fixture

`representatives.bam` contains three primary TAGET polished representative
alignments at SMaHT COLO829 SNV `chr1:948136:G>A` on GRCh38 no-alt:

| Representative | Base | Sequence allele | TAGET COLO829T count | TAGET COLO829BL count |
|---|---|---:|---:|---:|
| `PB.7.4` | G | 0 | 68 | 99 |
| `PB.7.5` | A | 1 | 2 | 0 |
| `PB.7.8` | A | 1 | 7 | 0 |

The sequences come from TAGET Zenodo record 8319497 `final.rep.fq` (SHA-256
`6b008d7610c03bdcc62e716bd365f4b448badc9c77554c58574068c7c9ab487b`).
The counts come from its `final.mapped_fl_count.txt` (SHA-256
`c54e284479ddb78714a30e03dfdec04e38b2d45fbb26f0a96717d83e5e9d80d9`).
The SNV comes from the SMaHT COLO829 SNV truth set (SHA-256
`e57f77c8bc280bbad2255856c4062348f92f25da38726f367322c52736578f31`).
The alignment used NCBI GRCh38 no-alt and minimap2 2.26
`-ax splice:hq -uf`. The fixture is a subset of the resulting sorted BAM;
its SHA-256 is `b94f58a456eb6b3533cd8803dc368435b7f9b3d924f418be9c15f4788c249868`.

These are polished representative sequences. The full-length counts describe
isoform clusters and are not original reads proven to cross the SNV.
