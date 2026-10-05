# Focused Sid vaccine and tumor RNA read subsets

This audit acquires regional original reads for **51 variants with a vaccine
membership claim, four additional RNA controls, and five fusion geometries** in
four pinned RNA products. The index-count, vaccine-overlap and source-variant
claims are unioned and retained independently. NR2F2 is both a vaccine claim and
a discordance control. These are historical metadata claims, not newly selected
treatments or demonstrated antigen presentation.

The fusion shortlist is TPST1::CRCP, GABBR1::SLC29A1,
PARD3B::CDKN2B-AS1/CDKN2B and OTUD7A::FMN1. The two distinct PARD3B geometries
remain separate, as do original nomination aliases. The RNA products are bulk T2
SARC0277, tagged ONT T2, Kamil 10x T2 and PacBio T1. Support is scoped to each
product and is never summed across these products or treated as independent
biological replicates.

The vaccine/control panels retain the following original records. These are
counts of acquired context across all 55 loci, not alternate-allele support:

| Product | Original records | Acquisition status | Timed-out linked queries |
| --- | ---: | --- | ---: |
| SARC0277 bulk RNA | 132,224 | incomplete | 3 |
| Tagged ONT T2 | 197,191 | incomplete | 70 |
| Kamil 10x T2 | 117,060 | bounded | 0 |
| PacBio T1 | 12,057 | incomplete | 1 |

The source metadata claims and discrepancies remain in the frozen inventory.
The 458,532 records above are a storage/input total across four products; they
are not pooled biological support. Original qualities and tags are retained.

The selected denominator is **220 small-variant/product pairs and 20
geometry/product pairs (40 fusion views)**. This is separate from the complete
catalogue's 173,866 nomination/product pairs. The full catalogue remains pending.
No unselected target becomes a negative RNA result.

Small-variant requests use the original VCF allele interval with 150-base
flanks. Fusion requests retain the existing breakend windows and pinned
transcript exons. Osteosarc recovers linked original mates and supplementary
records, preserving SAM payloads, qualities and labels. Oversized small-variant
requests partition without dropping targets. Budgets and failed recovery remain
explicit; acquired context counts are not mutation-support counts or VAFs.

The vaccine/control metadata snapshot is
`efb65d4b683bda162879c86c0ea889afcc64704b5228224c061109f306503b91`.
The SV inventory and Ensembl 115 references preserve the original audit pins.
The separately pinned dense ::IGKC validation scans all 773,288 original records
per orientation. Its candidate intervals differ from the earlier bounded run;
counts from different ORF intervals are not directly comparable.
The priority correction restores all 487 original records from the earlier
short-hypothesis witnesses to the retained provenance; the initial dense run
retained none. Independent SAM/translation checks verify a forward candidate in
279 fragments across 12 labels and a reverse candidate in 22 fragments across
three labels. These different candidate sequences do not establish translation.

## Reproduce and preserve

Use a durable audit directory, not `/tmp`. Copy the frozen inventory and
reference files into it, then run from the repository checkout:

```sh
python -m examples.acquire_sid_priority_reads AUDIT \
  --metadata-cache METADATA_CACHE --cache AUDIT/cache --workers 2
python -m examples.check_sid_priority_reads AUDIT OUTPUT
```

`priority-selection.json` records every allele, vaccine claim, selection reason,
geometry and product. Original BAMs, indexes, acquisition receipts and shared
provenance assets remain under the audit cache. Checkpoints bind input, code,
version and settings; changed runs require fresh result directories. The checker
independently recounts original BAM records, validates acquisition partitions and
source/input pins, and checks representative RNA witnesses and translations.
Keep the cache and pinned assets together when moving or backing up the subsets.
The published relative file pins can verify a relocated copy; resumable replay
still uses original receipt paths. Portable replay remains tracked in #440.

The current durable working directory is
`/Users/iskander/code/isovar/.cache/sid-priority-local-20261005`.
The October 5 run is checked in the [receipt ledger](data/sid-priority/priority-subsets.json.gz)
and [original BAM file index](data/sid-priority/input-files.tsv). All 220 selected
small-variant/product pairs and 40 fusion views are accounted for. Every
completed dense reconstruction scans exactly its independently recounted BAM
record total. The 197 pinned input/provenance/result files total **524 MiB**;
the working cache, including intermediate queries and preserved prior attempts,
occupies approximately **851 MB**. Original BAMs stay in the durable directory,
not Git. No whole library was downloaded.

The fusion run returns 193 event-related ORF hypotheses: 18 views have
event-linked candidates, four have splice-ambiguous candidates, and 18 have no
candidate paths under the recorded search settings. Five representatives, one
per geometry, pass independent original-query sequence, translation, fragment,
label, missing-quality and original-BAM witness checks. These checks establish
RNA sequence/count evidence; they do not establish translation or SV origin.
Upstream incomplete recovery and bounded discovery remain visible in every view.

The complete catalogue and the all-source vaccine allele/protein audit in #218
remain pending. This focused acquisition does not convert those unexamined pairs
to negative results. Earlier engine results and sandbox network-denied attempts
are preserved separately and excluded from the checked current ledger.

The [non-SNV additions](sid-neoorf-read-subsets.md) extend the regional inputs to
the remaining indel/splice nominations and FOXO3/ATP5MG RNA candidates. They keep
this original selection and archive unchanged; broader neoORF discovery and the
small-variant allele/protein audit remain separate work.
