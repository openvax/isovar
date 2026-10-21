# COLO829 full-run selected-site fixture

`SRR8615617.all_sites.grch38.bam` contains original read alignments retained
after the entire 115,269,780-spot (230,539,560-read) SRA run was decoded and
aligned to GRCh38. The filtering BED covers every SMaHT COLO829 SNV/indel
source locus and 500 bases around both benchmark SV ends. This 10 MB BAM is
indexed and sufficient to regenerate the atlas tables offline; it is not a
whole-genome RNA alignment and cannot be used to discover unlisted events.
The BAM contains 129,534 alignment records, of which 125,980 are primary.

The two VCFs are unmodified public SV benchmark files. Their exact source
links, methods, and SHA-256 values are in
[`reports/hcc1395_colo829/README.md`](../../../reports/hcc1395_colo829/README.md).

| File | SHA-256 |
|---|---|
| BAM | `cbe534d21f80b5ffab35a55c6edbd4cdddf80f1c8505ec9b8d2e4b36867b3f44` |
| BAM index | `fac75a006a0c7681ea39dffd05186930d308cc38849e3d1d32047c3bd95fd8f2` |
| SV v2 VCF | `8ad075bf43aa71d89c591396569e813959ada44a30f171982c2d299411290bac` |
| SV 2022 lifted VCF | `a8f38454367cf8837a7dba1d29c9bc338a8e62cadfee6d2489931f8105f0dcdc` |

See [`COLO829_ATLAS.md`](../../../reports/hcc1395_colo829/COLO829_ATLAS.md)
for the source accession, SRA archive checksum, tool versions, alignment
commands, evidence definitions, and culture limitations.
