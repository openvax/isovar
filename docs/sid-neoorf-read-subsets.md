# Additional Sid non-SNV and fusion RNA inputs

The original [focused archive](sid-priority-read-subsets.md) already includes
seven indel nominations: GLIS3 chr9:3856149, H1_2, MAP2 chr2:209694768, MYO15B,
PIP5K1A, RNF213 and TECPR1. Its small-variant panel is an input acquisition,
not a completed allele/protein audit. Its 193 event-related ORF hypotheses
come from five fusion geometries, not from all non-SNV loci.

The follow-up adds all remaining indel and splice-site nominations in the same
frozen metadata, plus the documented FOXO3–STRADA/CCDC47 and ATP5MG–KMT2A RNA
breakpoints. That is **25 additional small-variant nominations and two fusion
geometries** in the same four products: bulk SARC0277, tagged ONT T2, Kamil 10x
T2 and PacBio T1. Of these 25 alleles, 24 are ready GRCh38 inputs; MUC3A's
nonliteral `dup` allele remains explicitly unassessable in every product.

KTN1, GTF3C5, the second GLIS3 nomination, the second MAP2 nomination and
COL3A1-Splice are included. Distinct nominations and allele placements remain
separate; their supports must not be added as independent mutations. In-frame
indels are included as sequence-context controls, not labeled frameshift ORFs.
Splice-site metadata nominates a locus to inspect; it does not establish an
altered splice junction or translated ORF.

The additional denominator is **100 small-variant/product pairs**, including
four unassessable MUC3A pairs, and **eight fusion geometry/product pairs
(16 oriented views)**. Together with the original selection, the inventory
accounts for 80 small-variant nominations and seven fusion geometries in four
products. Processing products and timepoints remain separate.

ATP5MG–KMT2A is absent from the original frozen SV inventory. A derived inventory
retains that original inventory and pins the existing original-RNA fixture
`tests/data/fusions/coding-corpus/ATP5MG--KMT2A.input.json.gz` as a supplemental
nomination. Its interbase donor/acceptor coordinates determine the retained
sides. The exact FOXO3 geometry is already in the original catalogue and matches
the existing ONT RNA fixture. See [#446](https://github.com/openvax/isovar/issues/446).
Neither nomination nor an observed RNA join demonstrates DNA-SV origin;
ATP5MG–KMT2A is compatible with read-through splicing.

## Reproduce the additions

Use the original audit as immutable input and a new durable directory:

```sh
python -m examples.acquire_sid_neoorf_reads PRIOR_AUDIT ADDITIONS \
  --metadata-cache METADATA_CACHE --cache ADDITIONS/cache --workers 2
python -m examples.check_sid_priority_reads ADDITIONS ADDITIONS/report
```

The script requires the same metadata snapshot as the original selection. It
subtracts already selected alleles and geometries, pins the parent selection,
parent inventory and supplemental fixture, then builds Ensembl 115 references.
Requests use original allele intervals with 150-base flanks, or both breakend
windows and transcript exons. Original sequence, qualities and tags are retained
with bounded mate/supplementary recovery and the same 8 GiB free-space guard.
No whole library is downloaded. The checker validates native BAM record totals,
the complete additional denominator and representative fusion witnesses.

These regional inputs support subsequent allele, splice and protein audits.
General cryptic initiation, intron-retention and exonization discovery remains
pending, as does the full 1,493-geometry SV catalogue. Acquired context records
are not alternate-allele counts, pooled biological support or proof of antigen
presentation. The [existing non-SNV report](../figures/osteosarc/neo-orfs-and-long-reads.md)
describes the earlier selected-fixture evidence and its limitations.
