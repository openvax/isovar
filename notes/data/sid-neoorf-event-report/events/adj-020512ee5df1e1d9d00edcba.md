# TRDV1::TRAC / shared:SV0631

Event ID: `adj-020512ee5df1e1d9d00edcba`.

[Report and counting definitions](../../../sid-neoorf-event-report.md). [Complete JSON ledger](../events.json.gz) retains all placements and provenance.

## Original DNA event

Origin status: `unknown_RNA_only_nomination`.

| Side | Nominated interbase cut | Retained side | Adjacent retained base (one-based) |
| --- | --- | --- | --- |
| 1 | chr14:22096638 | left | 22096638 |
| 2 | chr14:22547505 | right | 22547506 |

No original DNA VCF call is linked to this selected nomination. The RNA join alone does not identify a somatic DNA rearrangement.

Nomination aliases: `shared:SV0631`, `fusion-call:0169`.

## Transcript anatomy and event location

Transcript coordinates are in 5′→3′ transcript order. cDNA intervals are zero-based, half-open. Introns below are genomic one-based inclusive. The model does not establish a complete altered mature transcript.

| Side | Gene / transcript | Strand | Placement |
| --- | --- | --- | --- |
| 2 | TRAC / ENST00000611116 | + | exon 1/4 exonic; CDS unresolved/noncoding; cDNA [0,1) |


## RNA screen outcomes

Side 1 has no overlapping supplied Ensembl 115 transcript at its adjacent retained base; gene/CDS assignment there is unresolved.

| Product | Orientation | Status | Records / paths / event paths | ORF hypotheses |
| --- | --- | --- | --- | --- |
| 10x Kamil T2 | forward | splice_ambiguous_candidates | 17077 / 33 / 33 | 13 |
| 10x Kamil T2 | reverse | splice_ambiguous_candidates | 17077 / 33 / 33 | 0 |

## Every protein sequence / coding window

Each sequence is an RNA-supported hypothesis or local coding window. A terminal `*` means an in-window stop, not a demonstrated full-length protein. SV initiation/frame/polarity may be unresolved. In-frame windows are controls. Counts below are never summed across products, starts, windows or geometry aliases.

### seq-78a9de1c24c52de270d4

SV_RNA_ORF; **14 aa**; stop: **True**.

```text
MKQISRSLTLPCTS*
```

| Product / orientation | Complete fragments | Q20 complete | Contributing allele fragments | Barcode labels | Missing qualities |
| --- | --- | --- | --- | --- | --- |
| 10x Kamil T2 / forward | 1 | not assessed | see junction/path ledger | 1 | 0 |

Distinct RNA geometric placements (equivalent path/query placements consolidated here; every occurrence remains in the ledger). Blocks describe this ORF interval, not a proven full mature transcript:

| ORF genomic blocks in RNA order | Crossed RNA join(s) | Frame / initiation |
| --- | --- | --- |
| chr14:22096406–22096414(+) → chr14:22547506–22547541(+) | chr14:22096414(+) → chr14:22547506(+); unplaced `∅`; splice_ambiguous_event_junction | unresolved / consistent; tier sequence_only; noncoding_transcript; initiation unobserved |

Flags: `['annotated_path_frame_unresolved_or_ambiguous', 'initiation_unobserved', 'interval_base_quality_unassessed', 'no_annotated_start_in_supplied_models', 'peptide_novelty_unassessed', 'reconstruction_limit:extension_segment_limit', 'reverse_complement_support_only', 'rna_strand_unresolved', 'signal_lineage_unresolved', 'single_full_fragment_witness', 'splice_ambiguous_event_junction', 'translation_unobserved']`.

### seq-0d0c9942431124392ce2

SV_RNA_ORF; **156 aa**; stop: **True**.

```text
MPVRKAVTLNCLYETNIQNPDPAVYQLRDSKSSDKSVCLFTDFDSQTNVSQSKDSDVYITDKTVLDMRSMDFKSNSAVAW
SNKSDFACANAFNNSIIPEDTFFPSPESSCDVKLVEKSFETDTNLNFQNLSVIGFRILLLKVAGFNLLMTLRLWSS*
```

| Product / orientation | Complete fragments | Q20 complete | Contributing allele fragments | Barcode labels | Missing qualities |
| --- | --- | --- | --- | --- | --- |
| 10x Kamil T2 / forward | 0 | not assessed | see junction/path ledger | 0 | 0 |

Distinct RNA geometric placements (equivalent path/query placements consolidated here; every occurrence remains in the ledger). Blocks describe this ORF interval, not a proven full mature transcript:

| ORF genomic blocks in RNA order | Crossed RNA join(s) | Frame / initiation |
| --- | --- | --- |
| chr14:22096369–22096414(+) → chr14:22547506–22547778(+) → chr14:22549638–22549682(+) → chr14:22550557–22550663(+) | chr14:22096414(+) → chr14:22547506(+); unplaced `∅`; splice_ambiguous_event_junction | unresolved / consistent; tier sequence_only; noncoding_transcript; initiation unobserved |

Flags: `['annotated_path_frame_unresolved_or_ambiguous', 'initiation_unobserved', 'interval_base_quality_unassessed', 'no_annotated_start_in_supplied_models', 'no_full_fragment_witness', 'peptide_novelty_unassessed', 'reconstruction_limit:extension_segment_limit', 'rna_strand_unresolved', 'splice_ambiguous_event_junction', 'translation_unobserved']`.

### seq-14a21b70eb4708944bc1

SV_RNA_ORF; **188 aa**; stop: **True**.

```text
MLFSSLLCVFVAFSYSGSSVAQKVTQAQSSVSMPVRKAVTLNCLYETNIQKPDPAVYQLRDSKSSDKSVCLFTDFDSQTN
VSQSKDSDVYITDKTVLDMRSMDFKSNSAVAWSNKSDFACANAFNNSIIPADTFFPSPESSCDVKLVEKSFETDTNLNFQ
NLSVIGFRILLLKVAGFNLLMTLRLWSS*
```

| Product / orientation | Complete fragments | Q20 complete | Contributing allele fragments | Barcode labels | Missing qualities |
| --- | --- | --- | --- | --- | --- |
| 10x Kamil T2 / forward | 0 | not assessed | see junction/path ledger | 0 | 0 |

Distinct RNA geometric placements (equivalent path/query placements consolidated here; every occurrence remains in the ledger). Blocks describe this ORF interval, not a proven full mature transcript:

| ORF genomic blocks in RNA order | Crossed RNA join(s) | Frame / initiation |
| --- | --- | --- |
| chr14:22096050–22096098(+) → chr14:22096322–22096414(+) → chr14:22547506–22547778(+) → chr14:22549638–22549682(+) → chr14:22550557–22550663(+) | chr14:22096414(+) → chr14:22547506(+); unplaced `∅`; splice_ambiguous_event_junction | unresolved / consistent; tier sequence_only; noncoding_transcript; initiation unobserved |

Flags: `['annotated_path_frame_unresolved_or_ambiguous', 'initiation_unobserved', 'interval_base_quality_unassessed', 'no_annotated_start_in_supplied_models', 'no_full_fragment_witness', 'peptide_novelty_unassessed', 'reconstruction_limit:extension_segment_limit', 'rna_strand_unresolved', 'splice_ambiguous_event_junction', 'translation_unobserved']`.

### seq-1b830f2e0cce8b144fa4

SV_RNA_ORF; **27 aa**; stop: **True**.

```text
MKPTLISVLVIIFILNIQNPDPAVYQL*
```

| Product / orientation | Complete fragments | Q20 complete | Contributing allele fragments | Barcode labels | Missing qualities |
| --- | --- | --- | --- | --- | --- |
| 10x Kamil T2 / forward | 0 | not assessed | see junction/path ledger | 0 | 0 |

Distinct RNA geometric placements (equivalent path/query placements consolidated here; every occurrence remains in the ledger). Blocks describe this ORF interval, not a proven full mature transcript:

| ORF genomic blocks in RNA order | Crossed RNA join(s) | Frame / initiation |
| --- | --- | --- |
| chr14:21990496–21990541(+) → chr14:22547506–22547543(+) | chr14:21990541(+) → chr14:22547506(+); unplaced `∅`; splice_ambiguous_event_junction | unresolved / consistent; tier sequence_only; outside_supplied_transcripts; initiation unobserved |

Flags: `['annotated_path_frame_unresolved_or_ambiguous', 'initiation_unobserved', 'interval_base_quality_unassessed', 'no_annotated_start_in_supplied_models', 'no_full_fragment_witness', 'peptide_novelty_unassessed', 'reconstruction_limit:extension_segment_limit', 'rna_strand_unresolved', 'splice_ambiguous_event_junction', 'translation_unobserved']`.

### seq-1fb1353650133d369bf1

SV_RNA_ORF; **188 aa**; stop: **True**.

```text
MLFSSLLCVFVAFSYSGSSVAQKVTQAQSSVSMPVRKAVTLNCLYETNIQNPDPAVYQLRDSKSSDKSVCLFTDFDSQTN
VSQSKDSDVYITDKTVLDMRSMDFKSNSAVAWSNKSDFACANAFNNSIIPEDTFFPSPESSCDVKLVEKSFETDTNLNFQ
NLSVIGFRILLLKVAGFNLLMTLRLWSS*
```

| Product / orientation | Complete fragments | Q20 complete | Contributing allele fragments | Barcode labels | Missing qualities |
| --- | --- | --- | --- | --- | --- |
| 10x Kamil T2 / forward | 0 | not assessed | see junction/path ledger | 0 | 0 |

Distinct RNA geometric placements (equivalent path/query placements consolidated here; every occurrence remains in the ledger). Blocks describe this ORF interval, not a proven full mature transcript:

| ORF genomic blocks in RNA order | Crossed RNA join(s) | Frame / initiation |
| --- | --- | --- |
| chr14:22096050–22096098(+) → chr14:22096322–22096414(+) → chr14:22547506–22547778(+) → chr14:22549638–22549682(+) → chr14:22550557–22550663(+) | chr14:22096414(+) → chr14:22547506(+); unplaced `∅`; splice_ambiguous_event_junction | unresolved / consistent; tier sequence_only; noncoding_transcript; initiation unobserved |

Flags: `['annotated_path_frame_unresolved_or_ambiguous', 'initiation_unobserved', 'interval_base_quality_unassessed', 'no_annotated_start_in_supplied_models', 'no_full_fragment_witness', 'peptide_novelty_unassessed', 'reconstruction_limit:extension_segment_limit', 'rna_strand_unresolved', 'splice_ambiguous_event_junction', 'translation_unobserved']`.

### seq-2f806ae48c624a4d4ba8

SV_RNA_ORF; **156 aa**; stop: **True**.

```text
MKPTLISVLVIIFILNIQNPDPAVYQLRDSKSSDKSVCLFTDFDSQTNVSQSKDSDVYITDKTVLDMRSMDFKSNSAVAW
SNKSDFACANAFNNSIIPADTFFPSPESSCDVKLVEKSFETDTNLNFQNLSVIGFRILLLKVAGFNLLMTLRLWSS*
```

| Product / orientation | Complete fragments | Q20 complete | Contributing allele fragments | Barcode labels | Missing qualities |
| --- | --- | --- | --- | --- | --- |
| 10x Kamil T2 / forward | 0 | not assessed | see junction/path ledger | 0 | 0 |

Distinct RNA geometric placements (equivalent path/query placements consolidated here; every occurrence remains in the ledger). Blocks describe this ORF interval, not a proven full mature transcript:

| ORF genomic blocks in RNA order | Crossed RNA join(s) | Frame / initiation |
| --- | --- | --- |
| chr14:21990496–21990541(+) → chr14:22547506–22547778(+) → chr14:22549638–22549682(+) → chr14:22550557–22550663(+) | chr14:21990541(+) → chr14:22547506(+); unplaced `∅`; splice_ambiguous_event_junction | unresolved / consistent; tier sequence_only; outside_supplied_transcripts; initiation unobserved |

Flags: `['annotated_path_frame_unresolved_or_ambiguous', 'initiation_unobserved', 'interval_base_quality_unassessed', 'no_annotated_start_in_supplied_models', 'no_full_fragment_witness', 'peptide_novelty_unassessed', 'reconstruction_limit:extension_segment_limit', 'rna_strand_unresolved', 'splice_ambiguous_event_junction', 'translation_unobserved']`.

### seq-35660af74ec4fafaca39

SV_RNA_ORF; **156 aa**; stop: **True**.

```text
MPVRKAVTLNCLYETNIQKPDPAVYQLRDSKSSDKSVCLFTDFDSQTNVSQSKDSDVYITDKTVLDMRSMDFKSNSAVAW
SNKSDFACANAFNNSIIPEDTFFPSPESSCDVKLVEKSFETDTNLNFQNLSVIGFRILLLKVAGFNLLMTLRLWSS*
```

| Product / orientation | Complete fragments | Q20 complete | Contributing allele fragments | Barcode labels | Missing qualities |
| --- | --- | --- | --- | --- | --- |
| 10x Kamil T2 / forward | 0 | not assessed | see junction/path ledger | 0 | 0 |

Distinct RNA geometric placements (equivalent path/query placements consolidated here; every occurrence remains in the ledger). Blocks describe this ORF interval, not a proven full mature transcript:

| ORF genomic blocks in RNA order | Crossed RNA join(s) | Frame / initiation |
| --- | --- | --- |
| chr14:22096369–22096414(+) → chr14:22547506–22547778(+) → chr14:22549638–22549682(+) → chr14:22550557–22550663(+) | chr14:22096414(+) → chr14:22547506(+); unplaced `∅`; splice_ambiguous_event_junction | unresolved / consistent; tier sequence_only; noncoding_transcript; initiation unobserved |

Flags: `['annotated_path_frame_unresolved_or_ambiguous', 'initiation_unobserved', 'interval_base_quality_unassessed', 'no_annotated_start_in_supplied_models', 'no_full_fragment_witness', 'peptide_novelty_unassessed', 'reconstruction_limit:extension_segment_limit', 'rna_strand_unresolved', 'splice_ambiguous_event_junction', 'translation_unobserved']`.

### seq-5e5aa093cf0da8d8b8eb

SV_RNA_ORF; **14 aa**; stop: **True**.

```text
MKQISRTLTLPCTS*
```

| Product / orientation | Complete fragments | Q20 complete | Contributing allele fragments | Barcode labels | Missing qualities |
| --- | --- | --- | --- | --- | --- |
| 10x Kamil T2 / forward | 0 | not assessed | see junction/path ledger | 0 | 0 |

Distinct RNA geometric placements (equivalent path/query placements consolidated here; every occurrence remains in the ledger). Blocks describe this ORF interval, not a proven full mature transcript:

| ORF genomic blocks in RNA order | Crossed RNA join(s) | Frame / initiation |
| --- | --- | --- |
| chr14:22096406–22096414(+) → chr14:22547506–22547541(+) | chr14:22096414(+) → chr14:22547506(+); unplaced `∅`; splice_ambiguous_event_junction | unresolved / consistent; tier sequence_only; noncoding_transcript; initiation unobserved |

Flags: `['annotated_path_frame_unresolved_or_ambiguous', 'initiation_unobserved', 'interval_base_quality_unassessed', 'no_annotated_start_in_supplied_models', 'no_full_fragment_witness', 'peptide_novelty_unassessed', 'reconstruction_limit:extension_segment_limit', 'rna_strand_unresolved', 'splice_ambiguous_event_junction', 'translation_unobserved']`.

### seq-61cf35d0c3c3a39f81a9

SV_RNA_ORF; **156 aa**; stop: **True**.

```text
MPVRKAVTLNCLYETNIQNPDPAVYQLRDSKSSDKSVCLFTDFDSQTNVSQSKDSDVYITDKTVLDMRSMDFKSNSAVAW
SNKSDFACANAFNNSIIPADTFFPSPESSCDVKLVEKSFETDTNLNFQNLSVIGFRILLLKVAGFNLLMTLRLWSS*
```

| Product / orientation | Complete fragments | Q20 complete | Contributing allele fragments | Barcode labels | Missing qualities |
| --- | --- | --- | --- | --- | --- |
| 10x Kamil T2 / forward | 0 | not assessed | see junction/path ledger | 0 | 0 |

Distinct RNA geometric placements (equivalent path/query placements consolidated here; every occurrence remains in the ledger). Blocks describe this ORF interval, not a proven full mature transcript:

| ORF genomic blocks in RNA order | Crossed RNA join(s) | Frame / initiation |
| --- | --- | --- |
| chr14:22096369–22096414(+) → chr14:22547506–22547778(+) → chr14:22549638–22549682(+) → chr14:22550557–22550663(+) | chr14:22096414(+) → chr14:22547506(+); unplaced `∅`; splice_ambiguous_event_junction | unresolved / consistent; tier sequence_only; noncoding_transcript; initiation unobserved |

Flags: `['annotated_path_frame_unresolved_or_ambiguous', 'initiation_unobserved', 'interval_base_quality_unassessed', 'no_annotated_start_in_supplied_models', 'no_full_fragment_witness', 'peptide_novelty_unassessed', 'reconstruction_limit:extension_segment_limit', 'rna_strand_unresolved', 'splice_ambiguous_event_junction', 'translation_unobserved']`.

### seq-a64aa4667b0340830348

SV_RNA_ORF; **156 aa**; stop: **True**.

```text
MKPTLISVLVIIFILNIQNPDPAVYQLRDSKSSDKSVCLFTDFDSQTNVSQSKDSDVYITDKTVLDMRSMDFKSNSAVAW
SNKSDFACANAFNNSIIPEDTFFPSPESSCDVKLVEKSFETDTNLNFQNLSVIGFRILLLKVAGFNLLMTLRLWSS*
```

| Product / orientation | Complete fragments | Q20 complete | Contributing allele fragments | Barcode labels | Missing qualities |
| --- | --- | --- | --- | --- | --- |
| 10x Kamil T2 / forward | 0 | not assessed | see junction/path ledger | 0 | 0 |

Distinct RNA geometric placements (equivalent path/query placements consolidated here; every occurrence remains in the ledger). Blocks describe this ORF interval, not a proven full mature transcript:

| ORF genomic blocks in RNA order | Crossed RNA join(s) | Frame / initiation |
| --- | --- | --- |
| chr14:21990496–21990541(+) → chr14:22547506–22547778(+) → chr14:22549638–22549682(+) → chr14:22550557–22550663(+) | chr14:21990541(+) → chr14:22547506(+); unplaced `∅`; splice_ambiguous_event_junction | unresolved / consistent; tier sequence_only; outside_supplied_transcripts; initiation unobserved |

Flags: `['annotated_path_frame_unresolved_or_ambiguous', 'initiation_unobserved', 'interval_base_quality_unassessed', 'no_annotated_start_in_supplied_models', 'no_full_fragment_witness', 'peptide_novelty_unassessed', 'reconstruction_limit:extension_segment_limit', 'rna_strand_unresolved', 'splice_ambiguous_event_junction', 'translation_unobserved']`.

### seq-a7f48f867385a91f169f

SV_RNA_ORF; **156 aa**; stop: **True**.

```text
MPVRKAVTLNCLYETNIQKPDPAVYQLRDSKSSDKSVCLFTDFDSQTNVSQSKDSDVYITDKTVLDMRSMDFKSNSAVAW
SNKSDFACANAFNNSIIPADTFFPSPESSCDVKLVEKSFETDTNLNFQNLSVIGFRILLLKVAGFNLLMTLRLWSS*
```

| Product / orientation | Complete fragments | Q20 complete | Contributing allele fragments | Barcode labels | Missing qualities |
| --- | --- | --- | --- | --- | --- |
| 10x Kamil T2 / forward | 0 | not assessed | see junction/path ledger | 0 | 0 |

Distinct RNA geometric placements (equivalent path/query placements consolidated here; every occurrence remains in the ledger). Blocks describe this ORF interval, not a proven full mature transcript:

| ORF genomic blocks in RNA order | Crossed RNA join(s) | Frame / initiation |
| --- | --- | --- |
| chr14:22096369–22096414(+) → chr14:22547506–22547778(+) → chr14:22549638–22549682(+) → chr14:22550557–22550663(+) | chr14:22096414(+) → chr14:22547506(+); unplaced `∅`; splice_ambiguous_event_junction | unresolved / consistent; tier sequence_only; noncoding_transcript; initiation unobserved |

Flags: `['annotated_path_frame_unresolved_or_ambiguous', 'initiation_unobserved', 'interval_base_quality_unassessed', 'no_annotated_start_in_supplied_models', 'no_full_fragment_witness', 'peptide_novelty_unassessed', 'reconstruction_limit:extension_segment_limit', 'rna_strand_unresolved', 'splice_ambiguous_event_junction', 'translation_unobserved']`.

### seq-b5819de2962b36629e71

SV_RNA_ORF; **188 aa**; stop: **True**.

```text
MLFSSLLCVFVAFSYSGSSVAQKVTQAQSSVSMPVRKAVTLNCLYETNIQKPDPAVYQLRDSKSSDKSVCLFTDFDSQTN
VSQSKDSDVYITDKTVLDMRSMDFKSNSAVAWSNKSDFACANAFNNSIIPEDTFFPSPESSCDVKLVEKSFETDTNLNFQ
NLSVIGFRILLLKVAGFNLLMTLRLWSS*
```

| Product / orientation | Complete fragments | Q20 complete | Contributing allele fragments | Barcode labels | Missing qualities |
| --- | --- | --- | --- | --- | --- |
| 10x Kamil T2 / forward | 0 | not assessed | see junction/path ledger | 0 | 0 |

Distinct RNA geometric placements (equivalent path/query placements consolidated here; every occurrence remains in the ledger). Blocks describe this ORF interval, not a proven full mature transcript:

| ORF genomic blocks in RNA order | Crossed RNA join(s) | Frame / initiation |
| --- | --- | --- |
| chr14:22096050–22096098(+) → chr14:22096322–22096414(+) → chr14:22547506–22547778(+) → chr14:22549638–22549682(+) → chr14:22550557–22550663(+) | chr14:22096414(+) → chr14:22547506(+); unplaced `∅`; splice_ambiguous_event_junction | unresolved / consistent; tier sequence_only; noncoding_transcript; initiation unobserved |

Flags: `['annotated_path_frame_unresolved_or_ambiguous', 'initiation_unobserved', 'interval_base_quality_unassessed', 'no_annotated_start_in_supplied_models', 'no_full_fragment_witness', 'peptide_novelty_unassessed', 'reconstruction_limit:extension_segment_limit', 'rna_strand_unresolved', 'splice_ambiguous_event_junction', 'translation_unobserved']`.

### seq-c508e21b2e9bf3b72f38

SV_RNA_ORF; **188 aa**; stop: **True**.

```text
MLFSSLLCVFVAFSYSGSSVAQKVTQAQSSVSMPVRKAVTLNCLYETNIQNPDPAVYQLRDSKSSDKSVCLFTDFDSQTN
VSQSKDSDVYITDKTVLDMRSMDFKSNSAVAWSNKSDFACANAFNNSIIPADTFFPSPESSCDVKLVEKSFETDTNLNFQ
NLSVIGFRILLLKVAGFNLLMTLRLWSS*
```

| Product / orientation | Complete fragments | Q20 complete | Contributing allele fragments | Barcode labels | Missing qualities |
| --- | --- | --- | --- | --- | --- |
| 10x Kamil T2 / forward | 0 | not assessed | see junction/path ledger | 0 | 0 |

Distinct RNA geometric placements (equivalent path/query placements consolidated here; every occurrence remains in the ledger). Blocks describe this ORF interval, not a proven full mature transcript:

| ORF genomic blocks in RNA order | Crossed RNA join(s) | Frame / initiation |
| --- | --- | --- |
| chr14:22096050–22096098(+) → chr14:22096322–22096414(+) → chr14:22547506–22547778(+) → chr14:22549638–22549682(+) → chr14:22550557–22550663(+) | chr14:22096414(+) → chr14:22547506(+); unplaced `∅`; splice_ambiguous_event_junction | unresolved / consistent; tier sequence_only; noncoding_transcript; initiation unobserved |

Flags: `['annotated_path_frame_unresolved_or_ambiguous', 'initiation_unobserved', 'interval_base_quality_unassessed', 'no_annotated_start_in_supplied_models', 'no_full_fragment_witness', 'peptide_novelty_unassessed', 'reconstruction_limit:extension_segment_limit', 'rna_strand_unresolved', 'splice_ambiguous_event_junction', 'translation_unobserved']`.
