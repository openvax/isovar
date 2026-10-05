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

The selected denominator is **220 small-variant/product pairs and 20
geometry/product pairs (40 fusion views)**. This is separate from the complete
catalogue's 173,866 nomination/product pairs. The full catalogue remains pending.
No unselected target becomes a negative RNA result.

Small-variant requests use the original VCF allele interval with 2,000-base
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

The current durable working directory is
`/Users/iskander/code/isovar/.cache/sid-priority-20261005`.
Acquisition and its final checked ledger are in progress. Measured size,
outcomes and independent-check totals will be published after verification.
