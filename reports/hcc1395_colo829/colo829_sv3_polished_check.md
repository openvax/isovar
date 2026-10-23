# COLO829 SV3: polished-isoform sequence cross-check

**No exact polished SV3 join or resolved proteoform.** Across
68,553 TAGET full-length representative sequences,
the forward 19-base junction motif occurs in 0 and
the five strand-aware 27-base junction windows have 0 matches.

12 representatives have a genomic
span across both nominated DNA breakend coordinates, but a genomic
span does not establish an exon or RNA junction at either end.

The reverse complement of the 19-base motif occurs in two
representatives. It also occurs in the ordinary GRCh38 sequence at
chr1:224598498. A reverse-motif hit alone is not an SV join.

| Representative | Tumor full-length count | Normal full-length count | Matched motif | Tumor GFF exons |
|---|---:|---:|---|---|
| PB.1153.13 | 3 | 0 | junction_reverse_complement | chr1:224529182-224529420; chr1:224529915-224530163; chr1:224566199-224566264; chr1:224583164-224583267; chr1:224597269-224598776 |
| PB.1153.23 | 4 | 0 | junction_reverse_complement | chr1:224529195-224529420; chr1:224546829-224546939; chr1:224566199-224566264; chr1:224583164-224583267; chr1:224597269-224598771 |

Both reverse-motif representatives have a final exon in the
downstream reference block and no exon at the nominated upstream
SV3 breakend. Their tumor cluster counts cannot be assigned to
the SV3 junction. The separate-culture polished data cannot
exclude a rare or unpolished SV transcript.

Sources: [TAGET Zenodo record](https://zenodo.org/records/10091914) and the
[source-pinned atlas](README.md#input-sha-256). Regenerate with
python -m reports.hcc1395_colo829.build_colo829_sv3_polished_check
--fastq /path/to/final.rep.fq --gff /path/to/final.COLO829T.gff.
