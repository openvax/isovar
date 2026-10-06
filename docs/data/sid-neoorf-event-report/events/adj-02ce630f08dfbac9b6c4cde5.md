# ::CNN2

Event ID: `adj-02ce630f08dfbac9b6c4cde5`.

[Report and counting definitions](../../../sid-neoorf-event-report.md). [Complete JSON ledger](../events.json.gz) retains all placements and provenance.

## Original DNA event

Origin status: `unknown_RNA_only_nomination`.

| Side | Nominated interbase cut | Retained side | Adjacent retained base (one-based) |
| --- | --- | --- | --- |
| 1 | chr19:1038882 | right | 1038883 |
| 2 | chr9:42971353 | right | 42971354 |

No original DNA VCF call is linked to this selected nomination. The RNA join alone does not identify a somatic DNA rearrangement.

Nomination aliases: `fusion-call:0605`, `fusion-call:1159`.

## Transcript anatomy and event location

Transcript coordinates are in 5′→3′ transcript order. cDNA intervals are zero-based, half-open. Introns below are genomic one-based inclusive. The model does not establish a complete altered mature transcript.

| Side | Gene / transcript | Strand | Placement |
| --- | --- | --- | --- |
| 1 | CNN2 / ENST00000263097 | + | exon 7/7 3′ UTR; cDNA [1966,1967) |
| 1 | CNN2 / ENST00000348419 | + | exon 6/6 3′ UTR; cDNA [1847,1848) |
| 1 | CNN2 / ENST00000562958 | + | exon 7/7 3′ UTR; cDNA [2008,2009) |
| 1 | CNN2 / ENST00000564572 | + | exon 1/1 exonic; CDS unresolved/noncoding; cDNA [2084,2085) |
| 1 | CNN2 / ENST00000565096 | + | exon 7/7 3′ UTR; cDNA [1923,1924) |
| 1 | CNN2 / ENST00000566695 | + | exon 6/6 3′ UTR; cDNA [1977,1978) |
| 1 | CNN2 / ENST00000865315 | + | exon 6/6 3′ UTR; cDNA [1841,1842) |
| 1 | CNN2 / ENST00000865316 | + | exon 5/5 3′ UTR; cDNA [1702,1703) |
| 1 | CNN2 / ENST00000926772 | + | exon 7/7 3′ UTR; cDNA [2021,2022) |
| 1 | CNN2 / ENST00000926773 | + | exon 4/4 3′ UTR; cDNA [1675,1676) |
| 1 | CNN2 / ENST00000926774 | + | exon 6/6 3′ UTR; cDNA [1835,1836) |
| 1 | CNN2 / ENST00000926775 | + | exon 6/6 3′ UTR; cDNA [1822,1823) |
| 1 | CNN2 / ENST00000926776 | + | exon 7/7 3′ UTR; cDNA [1945,1946) |
| 1 | CNN2 / ENST00000926777 | + | exon 7/7 3′ UTR; cDNA [1838,1839) |
| 1 | CNN2 / ENST00000926778 | + | exon 5/5 3′ UTR; cDNA [1695,1696) |
| 1 | CNN2 / ENST00000926779 | + | exon 7/7 3′ UTR; cDNA [1840,1841) |


## RNA screen outcomes

Side 2 has no overlapping supplied Ensembl 115 transcript at its adjacent retained base; gene/CDS assignment there is unresolved.

| Product | Orientation | Status | Records / paths / event paths | ORF hypotheses |
| --- | --- | --- | --- | --- |
| 10x Kamil T2 | forward | no_candidate_paths | 50000 / 0 / 0 | 0 |
| 10x Kamil T2 | reverse | no_candidate_paths | 50000 / 0 / 0 | 0 |

No resolved ORF/coding window in these acquired inputs and search bounds. For splice nominations, alternate exon use, intron retention and exonization have not been comprehensively reconstructed. Protein sequence and altered transcript structure remain unknown; reference anatomy above is the available description.
