# HCC1395 original-read SV junction fixture

This BAM contains 309 unmodified minimap2 alignment records from the
source-derived, 90,603-record selected-site alignment of the complete
HCC1395 ONT direct RNA run [SRR23456709](https://www.ncbi.nlm.nih.gov/sra/SRR23456709).
The input alignment is described in the
[whole-run fixture provenance](../hcc1395_whole_run/README.md). The BED selects
500 bp on either side of each nominated breakend for three candidate SVs.
The output is coordinate sorted and indexed. It is a regional fixture, not a
whole-genome RNA alignment or a representative estimate of expression.

Rebuild from the source-derived selected-site BAM:

```sh
samtools view -b -L tests/data/hcc1395_sv_junction/fixture_regions.bed \
  SRR23456709.all_sites.grch38.bam \
  -o tests/data/hcc1395_sv_junction/SRR23456709.sv_candidates.grch38.bam
samtools index tests/data/hcc1395_sv_junction/SRR23456709.sv_candidates.grch38.bam
```

SHA-256: BAM `d38efe0172a9d9e71c4ca47c89723d552b1956be98bdf8b2697d86bd41e187cd`;
BAI `0ecc584985738558413201979e225631f4b113c0edc1e7e1667d4da186ed9fa3`.

The three nominated adjacencies and their orientation are pinned in
`reports/hcc1395_colo829/build_sv_junction_audit.py`. Junction support counts
actual read segments that cross the join. One original read can support one
complete path, but the 3, 2 and 7 junction fragments are not 3, 2 and 7
independent transcript isoforms.

`ensembl87_models.json.gz` pins 27 GRCh38 transcript models (12 HNRNPK,
11 PKNOX1, 4 SUGT1) from Ensembl release 87. It is generated offline by
`python -m reports.hcc1395_colo829.build_sv_reference_models` using the
versioned Ensembl 87 cDNA and gene annotation already cached by PyEnsembl.
Every model has A/C/G/T cDNA and explicit 0-based half-open exons; complete
CDS coordinates are included only when the annotated frame validates.
No reference base is copied into an observed RNA path. SHA-256:
`8fc93e4457f993b1b6de9185512454ebfe5ac18b805dc060de372d7b6453f927`.
