# COLO829 SV3 local sequence fixtures

- reference_slice.json is the GRCh38 no-alt chr1 interval
  [224594494, 224598926), extracted from the NCBI FASTA with gzip SHA-256
  fb4243ebb014caf27111f24dd62b7ce42160f28581da6f8fcd6cba5977778d02.
- polished_reverse_motif_decoys.fasta contains 119-base windows centered on
  the reverse complement of the apparent SV3 unplaced sequence in TAGET
  representatives PB.1153.13 (polished offset 1890) and PB.1153.23
  (offset 1736). The original final.rep.fq has SHA-256
  6b008d7610c03bdcc62e716bd365f4b448badc9c77554c58574068c7c9ab487b
  and is available from [Zenodo 10091914](https://zenodo.org/records/10091914).
  These windows are negative controls: neither includes a plus-oriented
  exact SV3 junction 27-mer. The complete source is required for the
  68,553-representative cross-check.
