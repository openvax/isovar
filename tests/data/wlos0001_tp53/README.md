# WLOS_0001 canine TP53 source fixture

This 410 KB fixture pairs a **published tumor/normal DNA somatic call** with
**original aligned tumor RNA reads** from dog WLOS_0001. It is a single-locus
replay, not a whole-sample variant call set or a clinical interpretation.

Run from the repository root after installing Isovar's test dependencies:

```sh
python examples/wlos0001_tp53.py --output wlos0001-tp53-protein-hypotheses.json
```

The JSON and TSV retain every Isovar protein hypothesis and scoped RNA evidence.
The unit test checks allele and fragment support and a translated mutant Cys
hypothesis from original reads. RNA read counts are not independent molecules or
transcript abundance. The BAM retains PCR duplicate and secondary records;
Isovar's default read collector excludes them from the reported allele counts.

## Primary sources and sample match

- Gardner et al., *Communications Biology* 2:266 (2019), DOI
  [10.1038/s42003-019-0487-2](https://doi.org/10.1038/s42003-019-0487-2),
  Supplementary Data 1 maps dog 1 to `WLOS_0001` and lists tumor RNA,
  tumor DNA, and control DNA. Supplementary Data 6, `Primary Tumor` row for
  dog 1, reports `P53`, CanFam3.1 `5:32563056 G>A`, tumor DNA VAF 0.462,
  and `p.Arg263Cys` under that study's annotation. The local VCF records this
  DNA call; the normal DNA reads are not in this fixture.
- Original tumor RNA: [SRA SRR8691780](https://www.ncbi.nlm.nih.gov/sra/SRR8691780),
  BioProject [PRJNA525883](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA525883).
  Extracted with SRA Toolkit 3.4.1:

  ```sh
  sam-dump --aligned-region 5:32562800-32563300 --header SRR8691780
  ```

  The resulting 1,087 alignment records were coordinate-sorted and reheadered
  to retain chromosome 5 and the original read-group definitions. Read names,
  sequences, qualities, flags, tags, CIGAR strings, and original genomic
  coordinates were preserved. The BAM and BAI are local CI inputs.
- Reference: Ensembl release 98
  [CanFam3.1 chromosome 5 FASTA](https://ftp.ensembl.org/pub/release-98/fasta/canis_familiaris/dna/Canis_familiaris.CanFam3.1.dna.chromosome.5.fa.gz),
  [GTF](https://ftp.ensembl.org/pub/release-98/gtf/canis_familiaris/Canis_familiaris.CanFam3.1.98.gtf.gz),
  [cDNA FASTA](https://ftp.ensembl.org/pub/release-98/fasta/canis_familiaris/cdna/Canis_familiaris.CanFam3.1.cdna.all.fa.gz),
  and [peptide FASTA](https://ftp.ensembl.org/pub/release-98/fasta/canis_familiaris/pep/Canis_familiaris.CanFam3.1.pep.all.fa.gz).
  The GTF retains every feature for gene `ENSCAFG00000016714` (`TP53`),
  and the small cDNA/peptide FASTAs retain its three release-98 transcripts.
  Their sequences were independently checked to be byte-identical to the
  corresponding Ensembl release-98 cDNA and peptide sequences.

The compressed chromosome FASTA preserves all 88,915,250 chromosome-5
positions and the exact 13,512 bases of the TP53 gene (`32560598–32574109`);
all other bases are `N` masked. This keeps the original VCF/BAM coordinates
and source gene sequence without packaging a full chromosome. Do not use it
for reference validation or analysis outside TP53.

The paper's `p.Arg263Cys` label is annotation dependent. With the pinned
Ensembl release-98 TP53 transcripts, Varcode predicts `p.R261C` for
`ENSCAFT00000026465`, `p.R325C` for `ENSCAFT00000081814`, and `p.R238C` for
`ENSCAFT00000086831`. The genomic `G>A` call and the observed RNA allele are
the stable comparison keys here; the test does not assert a single protein
residue number across annotations.

## Source file SHA-256

| Source | SHA-256 |
|---|---|
| Ensembl 98 chromosome 5 FASTA gzip | `021df6e198011d5adad87678f064abb04fe3ff76cfbc20f0ad0be1caf4baaca1` |
| Ensembl 98 GTF gzip | `5d567adae8bcefef2e21b7d8d47bbb3beb6f0a9d4ae302c074867c3bebfdab2d` |
| Ensembl 98 cDNA FASTA gzip | `76070a4d11b550faad896c8500b001c45bb56ca5e31da85f3300f124f5406376` |
| Ensembl 98 peptide FASTA gzip | `60ce869e46ff855cee1543b5272b002c8598005c6a0d0ef9774e9a0f9db41515` |

The committed files can be checked with `shasum -a 256
tests/data/wlos0001_tp53/*`; the FASTA, GTF, VCF, and BAM are source-derived
subsets described above.
