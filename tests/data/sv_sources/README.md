# Original SV source fixture

`HCC1395.svscope.pcr.vcf` is the unmodified VCF inside `Extend Data File1.zip`
from [SVScope supplementary data, Zenodo 18688406](https://zenodo.org/records/18688406).
The dataset is licensed CC BY 4.0. Attribution: Tu, Zhang et al.,
[Genome Biology 27, 186 (2026)](https://doi.org/10.1186/s13059-026-04076-0).

Original filename: `HCC1395_minda_ensemble.Filter.PCR_Validated.vcf`.
SHA-256: `892222606e3951f20a8608b28ff32e6dc320257d573c9b832108156e9a361424`.
The 1,118 records include 32 explicit PCR-positive entries. The filename does
not establish validation of the other records. Coordinates, malformed endpoints,
consensus sequence fields and source support identifiers remain unchanged.

The COLO829 fixture is already pinned under
[`../colo829_whole_run`](../colo829_whole_run/README.md). Its four embedded END
breakend alleles are interpreted only under an explicitly selected profile for
that exact original checksum; the file is never rewritten.
