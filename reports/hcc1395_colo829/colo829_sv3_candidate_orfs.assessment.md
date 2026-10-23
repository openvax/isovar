# COLO829 SV3: RNA ORFs and proteoform status

**Proteoform call: not established.** The sequences below are exploratory
ATG-to-stop translations of an assembled RNA path, not measured protein.

The SMaHT/Paulin Paz DNA call is `COLO829_SOMATIC_SV3_DEL`. Three distinct
original tumor-RNA fragments support the same breakpoint-compatible join
with `CTATAGTAGAACTTCTTTT` between its placed flanks.
Their minimum sequenced base qualities over these 19 bases are SRR8615617.101713327: 34, SRR8615617.7313798: 38, SRR8615617.97504216: 23.

Isovar finds 2 alternative ATG reading frames on a 1,213-base
assembled RNA context.
The 533 flanking voting fragments are regional
coverage; they are not support for either complete isoform. The assembly
reaches `extension_segment_limit`.

| Candidate | Amino acids | Shared join fragments | Minimum exact 27-nt witnesses | Whole ORF witnesses | Hypothetical unrearranged genomic ORF | Shared amino-acid prefix | Join-crossing 9-mers |
|---|---:|---:|---:|---:|---:|---:|---|
| 1 | 37 | 3 | 3 | 0 | 35 aa | 25 aa | `QAYSRTSFG`, `AYSRTSFGV` |
| 2 | 162 | 3 | 3 | 0 | 47 aa | 20 aa | `VRPIVELLL`, `RPIVELLLE`, `PIVELLLEL` |

An exact 27-nt witness confirms the nucleotide window
under a hypothetical peptide. It cannot select a reading frame
or establish translation.

## Candidate amino-acid sequences

Candidate 1 (`sv_orf_dff1d97445d3b77f86345b53426ac22fd72056fbd5929aa4cda0ef9f3bfe4966`):

```text
MCLPHVSVSYFGNSHNISNLPCCQAYSRTSFGVASYG
```

Candidate 2 (`sv_orf_ddeae530c3924c3fe7da51b5e5d090f07bcdfb853ece5143c45781b8047ee30a`):

```text
MSLCHILVILTIFQICPVVRPIVELLLELLLMDEQRKWFLEMKSTPGEDAMNIVELTTKDLQYYINLVNKVVAGYERIDPNFKRSSTIGKMLSNSIACYREIFHERKSQLRWQTSLLSYFRKLTQPPQPSATTTLISQQPSTCRQDPPPAKRLRLAEGSDDC
```

The unrearranged comparisons translate *intronic genomic sequence* from
the same ATGs. They are controls for sequence divergence, not known
reference proteins. None of the listed 9-mers occurs in the pinned Ensembl 87 annotated protein sequences.
This is not proteome-wide novelty for unannotated ORFs.

Ensembl 87 places the interval within the CNIH3 gene span, but no
annotated exon overlaps the assembled region. The bundled ORFdata
catalogue (2026-10-11.1) has no overlapping
coding block; optional ORF catalogues were not assessed. The selected-site
BAM cannot establish a mature splice path from this join to a CNIH3 exon.
Neither candidate has an annotated start or observed initiation. The
lack of a whole-ORF witness limits phasing; no single read is required
to span the assembled sequence. Ribosome and peptide data were not
assessed. A predicted proteoform would require a supported mature RNA
path and frame; protein evidence would validate its expression.

The existing `isovar sv-rna --orf-output-prefix` export writes the
[protein FASTA](colo829_sv3_candidate_orfs.protein.fasta),
[nucleotide FASTA](colo829_sv3_candidate_orfs.nucleotide.fasta),
[ORF JSON](colo829_sv3_candidate_orfs.json) and
[ORF TSV](colo829_sv3_candidate_orfs.tsv). This assessment is
[machine readable](colo829_sv3_candidate_orfs.assessment.json).
The independent [TAGET polished-sequence cross-check](colo829_sv3_polished_check.md) tests whether a full-length
representative contains the exact junction in either orientation.
To resolve one hypothesis, obtain original RNA reads through a
path-discriminating exon join and compare normal/reference splicing;
then seek frame-specific ribosome footprints or a junction-specific
peptide. The RNA and DNA came from separate COLO829 cultures.
