# GABBR1::SLC29A1 / shared:SV0111 / shared:SV0636

Event ID: `adj-52ace853d3c63e505e9e5bba`.

[Report and counting definitions](../../../sid-neoorf-event-report.md). [Complete JSON ledger](../events.json.gz) retains all placements and provenance.

## Original DNA event

Origin status: `DNA_VCF_linked`.

| Side | Nominated interbase cut | Retained side | Adjacent retained base (one-based) |
| --- | --- | --- | --- |
| 1 | chr6:29612905 | right | 29612906 |
| 2 | chr6:44218907 | right | 44218908 |

| DNA sample / caller | Original VCF allele | VF / SF / DF | PURPLE AF / JCN | Call / cluster |
| --- | --- | --- | --- | --- |
| ucla_T1_organoid / ESVEE_PURPLE | chr6:29612906 `C` → `[chr6:44218908[C` | 8 / 8 / 0 | 0.158,0.171 / 1.02 | 17723 / INV / COMPLEX |
| ucla_T2 / ESVEE_PURPLE | chr6:29612906 `C` → `[chr6:44218908[C` | 10 / 9 / 1 | 0.121,0.162 / 0.903 | 15486 / INV / PAIR_OTHER |
| ucla_T1 / ESVEE_PURPLE | chr6:29612906 `C` → `[chr6:44218908[C` | 12 / 12 / 0 | 0.161,0.205 / 1.07 | 14988 / INV / COMPLEX |

Original VCF sources:

- [ucla_T1_organoid](https://sid-sijbrandij-osteosarc-dataset.s3.us-west-2.amazonaws.com/kamil/oncoanalyser/IPISRC044_T1_organoid_ucla/purple/IPISRC044_tumor_T1_organoid_ucla.purple.sv.vcf.gz)
- [ucla_T2](https://sid-sijbrandij-osteosarc-dataset.s3.us-west-2.amazonaws.com/kamil/oncoanalyser/IPISRC044_T2_ucla/purple/IPISRC044_tumor_T2_ucla.purple.sv.vcf.gz)
- [ucla_T1](https://sid-sijbrandij-osteosarc-dataset.s3.us-west-2.amazonaws.com/kamil/oncoanalyser/IPISRC044_T1_ucla/purple/IPISRC044_tumor_T1_ucla.purple.sv.vcf.gz)

Nomination aliases: `shared:SV0111`, `shared:SV0636`, `fusion-call:1281`, `fusion-call:1293`, `fusion-call:1298`, `fusion-call:1366`, `fusion-call:1443`, `fusion-call:1531`.

## Transcript anatomy and event location

Transcript coordinates are in 5′→3′ transcript order. cDNA intervals are zero-based, half-open. Introns below are genomic one-based inclusive. The model does not establish a complete altered mature transcript.

| Side | Gene / transcript | Strand | Placement |
| --- | --- | --- | --- |
| 1 | GABBR1 / ENST00000355973 | - | intron chr6:29612615–29613242 |
| 1 | GABBR1 / ENST00000377012 | - | intron chr6:29612615–29613242 |
| 1 | GABBR1 / ENST00000377016 | - | intron chr6:29612615–29613242 |
| 1 | GABBR1 / ENST00000377034 | - | intron chr6:29612615–29613242 |
| 1 | GABBR1 / ENST00000462632 | - | intron chr6:29612615–29613242 |
| 1 | GABBR1 / ENST00000472823 | - | intron chr6:29612615–29613242 |
| 1 | GABBR1 / ENST00000476670 | - | intron chr6:29612615–29613242 |
| 1 | GABBR1 / ENST00000491829 | - | intron chr6:29612615–29613242 |
| 1 | GABBR1 / ENST00000494877 | - | intron chr6:29612615–29613242 |
| 1 | GABBR1 / ENST00000706533 | - | intron chr6:29612615–29613242 |
| 1 | GABBR1 / ENST00000967932 | - | intron chr6:29612615–29613242 |


## RNA screen outcomes

Side 2 has no overlapping supplied Ensembl 115 transcript at its adjacent retained base; gene/CDS assignment there is unresolved.

| Product | Orientation | Status | Records / paths / event paths | ORF hypotheses |
| --- | --- | --- | --- | --- |
| Bulk SARC0277 T2 | forward | no_candidate_paths | 12838 / 0 / 0 | 0 |
| Bulk SARC0277 T2 | reverse | no_candidate_paths | 12838 / 0 / 0 | 0 |
| ONT tagged T2 | forward | event_linked_candidates | 23204 / 11 / 11 | 3 |
| ONT tagged T2 | reverse | event_linked_candidates | 23204 / 11 / 11 | 2 |
| 10x Kamil T2 | forward | no_candidate_paths | 50000 / 0 / 0 | 0 |
| 10x Kamil T2 | reverse | no_candidate_paths | 50000 / 0 / 0 | 0 |
| PacBio T1 | forward | event_linked_candidates | 489 / 2 / 2 | 1 |
| PacBio T1 | reverse | event_linked_candidates | 489 / 2 / 2 | 1 |

## Every protein sequence / coding window

Each sequence is an RNA-supported hypothesis or local coding window. A terminal `*` means an in-window stop, not a demonstrated full-length protein. SV initiation/frame/polarity may be unresolved. In-frame windows are controls. Counts below are never summed across products, starts, windows or geometry aliases.

### seq-a96461edfbffe0aaf2d8

SV_RNA_ORF; **67 aa**; stop: **True**.

```text
MRSNLVSCPAWNSTATTSSSSLKDPGSRRPSWTSLAKERSQEQAKRNLEFQSPTLSPPMKATLSKPS*
```

| Product / orientation | Complete fragments | Q20 complete | Contributing allele fragments | Barcode labels | Missing qualities |
| --- | --- | --- | --- | --- | --- |
| PacBio T1 / forward | 9 | not assessed | see junction/path ledger | 8 | 9 |
| PacBio T1 / reverse | 9 | not assessed | see junction/path ledger | 8 | 9 |

Distinct RNA geometric placements (equivalent path/query placements consolidated here; every occurrence remains in the ledger). Blocks describe this ORF interval, not a proven full mature transcript:

| ORF genomic blocks in RNA order | Crossed RNA join(s) | Frame / initiation |
| --- | --- | --- |
| chr6:44230858–44230881(-) → chr6:44230810–44230889(+) → chr6:44231364–44231456(+) | chr6:44230858(-) → chr6:44230810(+); unplaced `CCCGCCT`; event_compatible_junction | reference_protein_only / consistent; tier sequence_only; antisense; initiation unobserved |

Flags: `['annotated_path_frame_unresolved_or_ambiguous', 'cell_umi_labels_unresolved', 'initiation_unobserved', 'interval_base_quality_unassessed', 'library_scope_unresolved', 'missing_base_qualities', 'no_annotated_start_in_supplied_models', 'peptide_novelty_unassessed', 'rna_strand_unresolved', 'signal_lineage_unresolved', 'translation_unobserved', 'unplaced_junction_sequence']`.

### seq-b46b49b974f528eaf2d8

SV_RNA_ORF; **37 aa**; stop: **True**.

```text
MKRLVSSSRAWWRMPVIPAPTEAEAGESLESGRRRLQ*
```

| Product / orientation | Complete fragments | Q20 complete | Contributing allele fragments | Barcode labels | Missing qualities |
| --- | --- | --- | --- | --- | --- |
| ONT tagged T2 / forward | 3 | not assessed | see junction/path ledger | 3 | 0 |

Distinct RNA geometric placements (equivalent path/query placements consolidated here; every occurrence remains in the ledger). Blocks describe this ORF interval, not a proven full mature transcript:

| ORF genomic blocks in RNA order | Crossed RNA join(s) | Frame / initiation |
| --- | --- | --- |
| chr6:29612906–29612925(-) → chr6:44218908–44219001(+) | chr6:29612906(-) → chr6:44218908(+); unplaced `∅`; breakpoint_junction | unresolved / consistent; tier sequence_only; intronic; initiation unobserved |
| chr6:29612906–29612925(-) → chr6:44218908–44219001(+) | chr6:29612906(-) → chr6:44218908(+); unplaced `∅`; breakpoint_junction | ambiguous / consistent; tier sequence_only; intronic; initiation unobserved |

Flags: `['annotated_path_frame_unresolved_or_ambiguous', 'initiation_unobserved', 'interval_base_quality_unassessed', 'library_scope_unresolved', 'no_annotated_start_in_supplied_models', 'path_end_truncated', 'peptide_novelty_unassessed', 'reconstruction_limit:extension_segment_limit', 'reconstruction_limit:input_acquisition_incomplete', 'reconstruction_limit:repeated_genomic_position', 'rna_strand_unresolved', 'signal_lineage_unresolved', 'translation_unobserved']`.

### seq-0b1babc187468405933d

SV_RNA_ORF; **237 aa**; stop: **True**.

```text
MRSNLVSCPAWNSTATTSSSSLKDPGSRRPSWTSLAKERSQEQAKESGVSVSNSQPTNESHSIKAILKNISVLAFSVCFI
FTITIGMFPAVTVEVKSSIAGSSTWERYFIPVSCFLTFNIFDWLGRSLTAVFMWPGKDSRWLPSLVLARLVFVPLLLLCN
IKPRRYLTVVFEHDAWFIFFMAAFAFSNGYLASLCMCFGPKKVKPAEAETAGAIMAFFLCLGLALGAVFSFLFRAIV*
```

| Product / orientation | Complete fragments | Q20 complete | Contributing allele fragments | Barcode labels | Missing qualities |
| --- | --- | --- | --- | --- | --- |
| ONT tagged T2 / forward | 1 | not assessed | see junction/path ledger | 1 | 0 |
| ONT tagged T2 / reverse | 1 | not assessed | see junction/path ledger | 1 | 0 |

Distinct RNA geometric placements (equivalent path/query placements consolidated here; every occurrence remains in the ledger). Blocks describe this ORF interval, not a proven full mature transcript:

| ORF genomic blocks in RNA order | Crossed RNA join(s) | Frame / initiation |
| --- | --- | --- |
| chr6:44230858–44230881(-) → chr6:44230810–44230889(+) → chr6:44231364–44231385(+) → chr6:44231388–44231461(+) → chr6:44231998–44232106(+) → chr6:44232343–44232428(+) → chr6:44232807–44233006(+) → chr6:44233417–44233528(+) | chr6:44230858(-) → chr6:44230810(+); unplaced `CCCGCCT`; event_compatible_junction | reference_protein_only / consistent; tier sequence_only; antisense; initiation unobserved |

Flags: `['annotated_path_frame_unresolved_or_ambiguous', 'initiation_unobserved', 'interval_base_quality_unassessed', 'library_scope_unresolved', 'no_annotated_start_in_supplied_models', 'peptide_novelty_unassessed', 'reconstruction_limit:extension_segment_limit', 'reconstruction_limit:input_acquisition_incomplete', 'reconstruction_limit:repeated_genomic_position', 'rna_strand_unresolved', 'signal_lineage_unresolved', 'single_full_fragment_witness', 'translation_unobserved', 'unplaced_junction_sequence']`.

### seq-819c9165d4cd333bf94f

SV_RNA_ORF; **34 aa**; stop: **True**.

```text
MRSNLVSCPAWNSTATTSSSSLKDPGAGDQVGPH*
```

| Product / orientation | Complete fragments | Q20 complete | Contributing allele fragments | Barcode labels | Missing qualities |
| --- | --- | --- | --- | --- | --- |
| ONT tagged T2 / forward | 1 | not assessed | see junction/path ledger | 1 | 0 |
| ONT tagged T2 / reverse | 1 | not assessed | see junction/path ledger | 1 | 0 |

Distinct RNA geometric placements (equivalent path/query placements consolidated here; every occurrence remains in the ledger). Blocks describe this ORF interval, not a proven full mature transcript:

| ORF genomic blocks in RNA order | Crossed RNA join(s) | Frame / initiation |
| --- | --- | --- |
| chr6:44230858–44230881(-) → chr6:44230810–44230852(+) → chr6:44230854–44230884(+) | chr6:44230858(-) → chr6:44230810(+); unplaced `CCCGCCT`; event_compatible_junction | reference_protein_only / consistent; tier sequence_only; antisense; initiation unobserved |

Flags: `['annotated_path_frame_unresolved_or_ambiguous', 'initiation_unobserved', 'interval_base_quality_unassessed', 'library_scope_unresolved', 'no_annotated_start_in_supplied_models', 'peptide_novelty_unassessed', 'reconstruction_limit:extension_segment_limit', 'reconstruction_limit:input_acquisition_incomplete', 'reconstruction_limit:repeated_genomic_position', 'rna_strand_unresolved', 'signal_lineage_unresolved', 'single_full_fragment_witness', 'translation_unobserved', 'unplaced_junction_sequence']`.
