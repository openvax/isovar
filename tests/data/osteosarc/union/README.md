# Full Sid small-variant RNA panel

This offline regression uses the original, allele-balanced T2 short-read
members of [openvax-v1](https://github.com/iskandr/osteosarc/releases/tag/openvax-v1).
Read selection remains owned by osteosarc. The 345 older Isovar members are
unchanged; the package adds 184 upstream members, without renaming or reselecting
their records. The manifest pins the shared bundle hash, target definitions,
member identities, independent record observations, and pipeline outcomes.

## Complete target accounting

The pinned recipe has **187 literal small variants**: 179 current catalog
alleles and eight additional historical, reference-specific or library alleles.
It also records **two unresolved catalog entries**. Structural targets are
outside this regression.

| Scope | Targets | Treatment |
|---|---:|---|
| Literal GRCh38 alleles | 184 | Every member runs through allele collection and the public RNA-to-protein pipeline |
| NR2F2 GRCh37; MT_ND5 GRCh37 and hg19 | 3 | No matching T2 member; keep their native coordinates and explicit status |
| MUC3A; OTUD4 | 2 | No literal allele in the pinned snapshot; retain the upstream reason |

The 184 members contain 4,323 records counted across members (not unique
molecules). Twenty-two members are empty selections. Counts describe these
selected fixtures, not sample depth, VAF, expression, or calibration controls.
Missing members, uncallable reads, missing coding models, and absence of an
alternate allele are distinct conditions. An exception always fails the test.

## Independent checks

`tests/osteosarc_union_helpers.py` reads the original CIGAR operations and
literal alleles. It checks reference/alternate/other evidence and query offsets,
including MNVs and complex replacements where both flanks are aligned. An `N`
operation is a reference skip, never evidence for a deletion. These direct
observations do not perform repeat normalization. The SNV tests additionally
compare the entire primary-record partition with `ReadCollector`, with mate
merging disabled and supplementary records explicitly excluded. The full
pipeline tests retain production defaults, including mate merging.

The reference subset is copied from original Ensembl GRCh38 release 87 GTF,
cDNA and peptide archives. `reference/manifest.json` pins their URLs and SHA-256
hashes and every subset asset. The existing independent reference builder maps
genomic coordinates through raw exon rows, applies the literal edit to cDNA,
and translates with NCBI genetic codes. It records incomplete/noncoding models
and transcript-specific exclusions rather than inventing coding sequences.
There are 571 validated models, providing coding expectations for 171 targets.
All 184 targets run, using the validated transcript whitelist where available.
This does not test alternative splice consequences for edits outside those
coding models.

Protein checks inspect every ranked candidate before the default top-one cap.
They independently verify reading frame, translated RNA sequence,
stop status, mutation interval and the nominated reference-edit peptide.
The frozen outcomes are:

| Pipeline outcome | Targets |
|---|---:|
| Protein matches the nominated reference edit | 34 |
| Protein contains additional RNA sequence differences | 2 |
| No alternate reads | 83 |
| No callable allele | 52 |
| Alternate reads, but no assembled RNA candidate | 10 |
| Alternate reads, but no validated coding model | 2 |
| RNA candidate, but no translation | 1 |

FCGBP and NBPF1 are the two additional-difference cases. Their observed RNA
translations and mutation metadata pass the independent checks; the exact
peptide differences from the nominated edit are retained in the ledger.
ANKRD20A8P and PCDH19 lack a validated coding expectation in this reference
subset. PPP1R3F has an RNA candidate without a returned translation. These are
recorded outcomes, not claims of absent expression or biological truth.

## Reproduce

From the repository root, after installing development dependencies:

```sh
python -m pytest tests/test_osteosarc_union.py tests/test_sid_data.py
python -m isovar.sid_data check --from isovar/data/sid-reads
python -m tests.data.osteosarc.union.build
```

The tests block socket connections, export from the packaged bundle into a
temporary cache, and index only the pinned local reference files. The builder
regenerates the ledger without acquiring any reads. To recreate the reference
subset, move the generated `reference` directory aside first, then run the
builder with `--reference-source /path/to/ensembl87`. It uses the three original
archives named in the reference manifest. Review regenerated expectations;
never accept a failing pipeline or translation check as a new baseline.

Scientific conventions: [SAM specification](https://samtools.github.io/hts-specs/SAMv1.pdf),
[Ensembl release 87](https://ftp.ensembl.org/pub/release-87/), and
[NCBI genetic codes](https://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi).
