# Mutation-supporting RNA in catalogue ncORF frames

Isovar 1.61.0 replays all seven ALT-supported ncORF-overlapping calls in the
COLO829/HCC1395 atlas. Nine supplied catalogue models bind to Ensembl 87 by
validating their entire ORF splice paths and proteins. Ordinary annotated CDS
models on those parent transcripts remain competitors in the export.

| Sample | ALT calls | ncORF models | Changed models, default coverage | Additional single-read hypotheses |
|---|---:|---:|---:|---:|
| COLO829 / SRR8615617 | 2 | 3 | 1 | 0 |
| HCC1395 / SRR23456709 | 5 | 6 | 0 | 2 |

COLO829 DYNLL2 `chr17:58089866:G>A` has 45 original ALT fragments. In the
`c17norep133` downstream ORF, the displayed 15-aa RNA sequence
`MPSLEALRCVCLWCV` differs from `MPSLGALRCVCLWCV` at G5E. It has 30 assembled
fragments, 21 fragments with one read spanning the entire coding window, and
minimum base depth 24. The reference ORF is 85 aa; the RNA candidate remains
partial. A longer 23-aa candidate has 28 assembled fragments, 11 single-read
spans and minimum base depth 15. The alternate DYNLL2 model and PITPNM1 model
are synonymous for their focal alleles and are not counted as mutant proteins.

HCC1395 MIR762HG and PPIP5K2 each yield a candidate supported by one original
ONT molecule only when minimum read coverage is reduced from two to one.
Nonfocal differences remain unclassified and may include sequencing errors.
NELFB affects the reference start codon; alternative initiation is unresolved.
UQCC1 and ST3GAL1 have no retained RNA/frame match under either budget. The
other MIR762HG isoform path is incompatible with the recovered reads.

## Sources and reproduction

`manifest.json` pins source BAMs, catalogue metadata, budgets and result hash.
`inputs/HCC1395.source.json` pins the original ALT-witness selection, SRA
recovery command, genome reference and all derived files. The 36 original
HCC1395 reads reproduce all five published ALT counts exactly. They are
selected witnesses from the same SRR23456709 library, not new replicates or an
allele-frequency denominator. Complete-genome alignment includes decoy/HLA
sequences; it is distinct from the older atlas no-alt reference.

With the pinned inputs, ORFdata 0.39.1 and cached Ensembl 87:

```sh
python -m reports.hcc1395_colo829.build_ncorf_rna_audit
python -m reports.hcc1395_colo829.build_evidence_pdfs
```

Rendering uses the committed `results.json.gz` and needs no remote RNA fetch.
RNA alternatives, model outcomes, sequences and scoped original witness IDs
remain in the result. Hundreds of overlapping local assembly windows are not
hundreds of independent ORFs. The two coverage passes reuse the same molecules.

The [library API](../../../docs/variant-rna-orfs.md) implements reconstruction;
report scripts select and display its results. Reference Ribo-seq/MS annotations
do not measure mutant initiation, translation or HLA presentation. Broader
intronic/cryptic-exon discovery remains tracked in
[#583](https://github.com/openvax/isovar/issues/583).
