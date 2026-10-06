# COL3A1

Event ID: `COL3A1-Splice`.

[Report and counting definitions](../../../sid-neoorf-event-report.md). [Complete JSON ledger](../events.json.gz) retains all placements and provenance.

## Original DNA event

Origin status: `catalogued_allele; original_DNA_support_unavailable`.

Original GRCh38 VCF-style allele: **chr2:189010889** `CGGTAGGAAACATTTTTCTCAATATAGGTCATAAAGCAGTCAGCATTTTAGTTTAATCATGCAAATTATTTTGAATAGAATAAATAAAATTAATAGGATGAAATAAAGATAGCATTTGGTATGAATTACATACATTGCATCTACTGATTCATTGCAGGGAATGTTAAATGCAACAAAATGGATCTTAGCCTCCAGATGAAAACCAGTTGAAATAAAAGCATTTAAAATTTAAATCCATAATGCAAACTTTATCAGATAATTGGGATAGTTACTATATGTTTTGAATAATACCTGTCACATTAACTCAGTTTGAAATATCTGTTTTTTAAAATTTTAAAGGTTAACTTCAAATCTCTCATTTGTTATTGTTATTTGTGAACAAAGGGAAAACTCACTTCTAATATTAGCGAATTTCATTGTGAGAGACCTATCCTCTTTTTAATAAACCATTTATAAACCTGTTAGCAATGGCCGGAACCAGGCCTCCTGAGGATGCACTGGTCTATAGCAATTGCCCTGCTGCATTTACTAAGAATCCTTACAAGGCATTTGTTTAAGAATATTGTTTATCAACTAAGAAGATTACAGCTTTGAAGTAGAGCAGGTCTCATATACATGAATAATAACATGGCACGATGAATGCTTCTTTAGAGTAAAAAGGTTTTCTTTAACTTGTTAAGTCAGAGTTGTCTAAGTAATTGTAATGTCATGATCATGTACATTTTGTCCTTTTTTACA` → `C`.

Source consequence: `unavailable`; source protein label: `unavailable`; source transcript: `NM_000090`; cDNA label: `c.4254+1_4255-1del`.

After removing VCF padding: genomic **[189010889,189011626)**; `GGTAGGAAACATTTTTCTCAATATAGGTCATAAAGCAGTCAGCATTTTAGTTTAATCATGCAAATTATTTTGAATAGAATAAATAAAATTAATAGGATGAAATAAAGATAGCATTTGGTATGAATTACATACATTGCATCTACTGATTCATTGCAGGGAATGTTAAATGCAACAAAATGGATCTTAGCCTCCAGATGAAAACCAGTTGAAATAAAAGCATTTAAAATTTAAATCCATAATGCAAACTTTATCAGATAATTGGGATAGTTACTATATGTTTTGAATAATACCTGTCACATTAACTCAGTTTGAAATATCTGTTTTTTAAAATTTTAAAGGTTAACTTCAAATCTCTCATTTGTTATTGTTATTTGTGAACAAAGGGAAAACTCACTTCTAATATTAGCGAATTTCATTGTGAGAGACCTATCCTCTTTTTAATAAACCATTTATAAACCTGTTAGCAATGGCCGGAACCAGGCCTCCTGAGGATGCACTGGTCTATAGCAATTGCCCTGCTGCATTTACTAAGAATCCTTACAAGGCATTTGTTTAAGAATATTGTTTATCAACTAAGAAGATTACAGCTTTGAAGTAGAGCAGGTCTCATATACATGAATAATAACATGGCACGATGAATGCTTCTTTAGAGTAAAAAGGTTTTCTTTAACTTGTTAAGTCAGAGTTGTCTAAGTAATTGTAATGTCATGATCATGTACATTTTGTCCTTTTTTACA` → `∅`; net length change **-737 nt**. This is genomic placement, not HGVS right normalization or a demonstrated mature-RNA edit.

Original DNA depth/VAF is unavailable in the frozen nomination.

## Transcript anatomy and event location

Transcript coordinates are in 5′→3′ transcript order. cDNA intervals are zero-based, half-open. Introns below are genomic one-based inclusive. The model does not establish a complete altered mature transcript.

| Side | Gene / transcript | Strand | Placement |
| --- | --- | --- | --- |
| allele | COL3A1 / ENST00000304636 | + | exon 50/51 CDS; cDNA [4370,4371); intron chr2:189010891–189011627 |
| allele | COL3A1 / ENST00000450867 | + | exon 49/50 CDS; cDNA [4271,4272); intron chr2:189010891–189011627 |
| allele | COL3A1 / ENST00000467886 | + | outside_transcript |
| allele | COL3A1 / ENST00000470167 | + | outside_transcript |
| allele | COL3A1 / ENST00000487010 | + | exon 2/3 exonic; CDS unresolved/noncoding; cDNA [1632,1633); intron chr2:189010891–189011627 |
| allele | COL3A1 / ENST00000713744 | + | exon 48/49 CDS; cDNA [4208,4209); intron chr2:189010891–189011627 |
| allele | COL3A1 / ENST00000713745 | + | exon 48/49 CDS; cDNA [4217,4218); intron chr2:189010891–189011627 |
| allele | COL3A1 / ENST00000879194 | + | exon 49/50 CDS; cDNA [4319,4320); intron chr2:189010891–189011627 |
| allele | COL3A1 / ENST00000879195 | + | exon 49/50 CDS; cDNA [4271,4272); intron chr2:189010891–189011627 |
| allele | COL3A1 / ENST00000879196 | + | exon 49/50 CDS; cDNA [4262,4263); intron chr2:189010891–189011627 |
| allele | COL3A1 / ENST00000879197 | + | exon 49/50 CDS; cDNA [4231,4232); intron chr2:189010891–189011627 |
| allele | COL3A1 / ENST00000879198 | + | exon 49/50 CDS; cDNA [4325,4326); intron chr2:189010891–189011627 |
| allele | COL3A1 / ENST00000879199 | + | exon 48/49 CDS; cDNA [4226,4227); intron chr2:189010891–189011627 |
| allele | COL3A1 / ENST00000879200 | + | exon 49/50 CDS; cDNA [4316,4317); intron chr2:189010891–189011627 |
| allele | COL3A1 / ENST00000879201 | + | exon 50/51 CDS; cDNA [4361,4362); intron chr2:189010891–189011627 |
| allele | COL3A1 / ENST00000879202 | + | exon 49/50 CDS; cDNA [4316,4317); intron chr2:189010891–189011627 |
| allele | COL3A1 / ENST00000879203 | + | exon 49/50 CDS; cDNA [4247,4248); intron chr2:189010891–189011627 |
| allele | COL3A1 / ENST00000957916 | + | exon 49/50 CDS; cDNA [4262,4263); intron chr2:189010891–189011627 |
| allele | COL3A1 / ENST00000957917 | + | exon 50/51 CDS; cDNA [4252,4253); intron chr2:189010891–189011627 |
| allele | COL3A1 / ENST00000957918 | + | exon 49/50 CDS; cDNA [4316,4317); intron chr2:189010891–189011627 |


## RNA screen outcomes

The event has no wholly exonic, non-boundary placement in the supplied models; the mature RNA edit is unknown. No resolved coding window establishes an altered protein in these inputs. An intronic or exon-boundary event can alter splicing; its mature-RNA outcome remains unresolved without an observed altered splice path. The genomic length change applies to a coding frame only where that literal edit is retained in CDS; intronic bases are not automatically translated.

| Product | Status | RNA ref / alt / other fragments | Coding windows |
| --- | --- | --- | --- |
| Bulk SARC0277 T2 | screened | 0 / 0 / 0 | 0 |
| ONT tagged T2 | screened | 0 / 0 / 25 | 0 |
| 10x Kamil T2 | screened | 0 / 0 / 0 | 0 |
| PacBio T1 | screened | 0 / 0 / 4 | 0 |

No resolved ORF/coding window in these acquired inputs and search bounds. For splice nominations, alternate exon use, intron retention and exonization have not been comprehensively reconstructed. Protein sequence and altered transcript structure remain unknown; reference anatomy above is the available description.
