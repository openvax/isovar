# HCC1395 whole-run original RNA regression fixture

`SRR23456709.fixture.grch38.bam` is a bounded selection of original ONT
direct RNA alignments from the **complete** 1,049,452-read SRA run
[SRR23456709](https://www.ncbi.nlm.nih.gov/sra/SRR23456709). The reads were
aligned to the NCBI GRCh38 no-alt analysis set with minimap2 2.26
`-ax splice -uf -k14`; `fixture_regions.bed` then selected genomic alignments
for OSTC, PRDX5, one translated indel, one reference-only call, one
no-overlap call, and two nominated SVs. The complete run was aligned before
any site selection. The BED and BAM are source-derived; they do not contain
manufactured allele sequences.

The validated original SRA archive SHA-256 is
`d69db8e48bbb20665c3e03a1fa35742be818a0127994e015cb3a560afae39ac8`.
The full selected-site BAM SHA-256 is
`9634b78414b9ddbbb9e098a01d899028dd0b2d703c3c9d0a5c24ce6a95cc9db9`.
This bounded BAM SHA-256 is
`c04f29e2e9f6b6fe719ae215ebb10d753aa724741c47ae0bf6f7b4a3c9978412`.
All BAMs are coordinate sorted, indexed, and pass `samtools quickcheck`.

Recreate the bounded BAM from the complete selected-site BAM:

```sh
samtools view -b -L tests/data/hcc1395_whole_run/fixture_regions.bed \
  SRR23456709.all_sites.grch38.bam \
  | samtools sort -o tests/data/hcc1395_whole_run/SRR23456709.fixture.grch38.bam -
samtools index tests/data/hcc1395_whole_run/SRR23456709.fixture.grch38.bam
```

The fixture reproduces read support and candidate-screen outcomes at its
selected loci. It cannot regenerate the complete 41,072-call atlas.
