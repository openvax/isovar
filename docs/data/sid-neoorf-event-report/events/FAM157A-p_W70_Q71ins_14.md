# FAM157A

Event ID: `FAM157A-p_W70_Q71ins_14`.

[Report and counting definitions](../../../sid-neoorf-event-report.md). [Complete JSON ledger](../events.json.gz) retains all placements and provenance.

## Original DNA event

Origin status: `catalogued_allele; original_DNA_support_unavailable`.

Original GRCh38 VCF-style allele: **chr3:198153259** `G` → `GGCGGCGGCGGCGGCAGCAGCAGCAGCAGCAGCAGCAGCAGCA`.

Source consequence: `unavailable`; source protein label: `unavailable`; source transcript: `NM_001145248`; cDNA label: `c.211_212ins(42)`.

After removing VCF padding: genomic **[198153259,198153259)**; `∅` → `GCGGCGGCGGCGGCAGCAGCAGCAGCAGCAGCAGCAGCAGCA`; net length change **+42 nt**. This is genomic placement, not HGVS right normalization or a demonstrated mature-RNA edit.

Original DNA depth/VAF is unavailable in the frozen nomination.

## Transcript anatomy and event location

Transcript coordinates are in 5′→3′ transcript order. cDNA intervals are zero-based, half-open. Introns below are genomic one-based inclusive. The model does not establish a complete altered mature transcript.

| Side | Gene / transcript | Strand | Placement |
| --- | --- | --- | --- |
| allele | ENSG00000290937 / ENST00000437428 | + | outside_transcript |
| allele | ENSG00000290937 / ENST00000634862 | + | outside_transcript |
| allele | ENSG00000290937 / ENST00000779226 | + | outside_transcript |
| allele | ENSG00000290937 / ENST00000779227 | + | exon 1/8 exonic; CDS unresolved/noncoding; cDNA [0,0); exon boundary |
| allele | ENSG00000290937 / ENST00000779228 | + | outside_transcript |
| allele | ENSG00000290937 / ENST00000779229 | + | exon 1/7 exonic; CDS unresolved/noncoding; cDNA [20,20) |
| allele | ENSG00000290937 / ENST00000779230 | + | outside_transcript |
| allele | ENSG00000290937 / ENST00000779231 | + | outside_transcript |
| allele | ENSG00000290937 / ENST00000779232 | + | outside_transcript |
| allele | ENSG00000290937 / ENST00000779233 | + | outside_transcript |
| allele | ENSG00000290937 / ENST00000779234 | + | outside_transcript |
| allele | ENSG00000290937 / ENST00000779235 | + | outside_transcript |
| allele | ENSG00000290937 / ENST00000779236 | + | outside_transcript |
| allele | ENSG00000290937 / ENST00000779237 | + | intron chr3:198152159–198167675 |
| allele | ENSG00000290937 / ENST00000779238 | + | outside_transcript |
| allele | ENSG00000290937 / ENST00000779239 | + | outside_transcript |
| allele | ENSG00000290937 / ENST00000779240 | + | outside_transcript |
| allele | ENSG00000290937 / ENST00000779241 | + | outside_transcript |
| allele | ENSG00000290937 / ENST00000779242 | + | outside_transcript |
| allele | ENSG00000290937 / ENST00000779243 | + | outside_transcript |
| allele | ENSG00000290937 / ENST00000779244 | + | intron chr3:198152159–198153323 |
| allele | ENSG00000290937 / ENST00000779245 | + | outside_transcript |
| allele | ENSG00000290937 / ENST00000779246 | + | exon 1/3 exonic; CDS unresolved/noncoding; cDNA [18,18) |
| allele | ENSG00000290937 / ENST00000779247 | + | outside_transcript |
| allele | ENSG00000290937 / ENST00000779248 | + | outside_transcript |
| allele | ENSG00000290937 / ENST00000779249 | + | outside_transcript |
| allele | ENSG00000290937 / ENST00000779250 | + | outside_transcript |
| allele | ENSG00000290937 / ENST00000779251 | + | outside_transcript |
| allele | ENSG00000290937 / ENST00000779252 | + | exon 2/4 exonic; CDS unresolved/noncoding; cDNA [256,256) |
| allele | ENSG00000290937 / ENST00000779253 | + | intron chr3:198152159–198153323 |
| allele | ENSG00000290937 / ENST00000779254 | + | intron chr3:198152159–198153323 |
| allele | ENSG00000290937 / ENST00000779255 | + | exon 2/4 exonic; CDS unresolved/noncoding; cDNA [199,199) |
| allele | ENSG00000290937 / ENST00000779256 | + | exon 2/3 exonic; CDS unresolved/noncoding; cDNA [195,195) |
| allele | ENSG00000290937 / ENST00000779257 | + | intron chr3:198152159–198153323 |
| allele | ENSG00000290937 / ENST00000779258 | + | exon 2/3 exonic; CDS unresolved/noncoding; cDNA [167,167) |
| allele | ENSG00000290937 / ENST00000779259 | + | exon 2/4 exonic; CDS unresolved/noncoding; cDNA [350,350) |
| allele | ENSG00000290937 / ENST00000779260 | + | exon 1/3 exonic; CDS unresolved/noncoding; cDNA [18,18) |
| allele | ENSG00000290937 / ENST00000779261 | + | exon 1/3 exonic; CDS unresolved/noncoding; cDNA [18,18) |
| allele | ENSG00000290937 / ENST00000779262 | + | outside_transcript |
| allele | ENSG00000290937 / ENST00000779263 | + | outside_transcript |
| allele | ENSG00000290937 / ENST00000779264 | + | outside_transcript |
| allele | ENSG00000301525 / ENST00000779477 | - | intron chr3:198136569–198160198 |
| allele | ENSG00000301525 / ENST00000779478 | - | outside_transcript |
| allele | ENSG00000301525 / ENST00000779479 | - | outside_transcript |
| allele | ENSG00000301525 / ENST00000779480 | - | outside_transcript |
| allele | ENSG00000301566 / ENST00000779731 | - | exon 1/2 exonic; CDS unresolved/noncoding; cDNA [111,111) |


## RNA screen outcomes

Plus-strand local allele, if the affected exon is retained: `∅` → `GCGGCGGCGGCGGCAGCAGCAGCAGCAGCAGCAGCAGCAGCA`; Minus-strand local allele, if the affected exon is retained: `∅` → `TGCTGCTGCTGCTGCTGCTGCTGCTGCTGCCGCCGCCGCCGC`. No resolved coding window establishes an altered protein in these inputs. An intronic or exon-boundary event can alter splicing; its mature-RNA outcome remains unresolved without an observed altered splice path. The genomic length change applies to a coding frame only where that literal edit is retained in CDS; intronic bases are not automatically translated.

| Product | Status | RNA ref / alt / other fragments | Coding windows |
| --- | --- | --- | --- |
| Bulk SARC0277 T2 | screened | 0 / 0 / 0 | 0 |
| ONT tagged T2 | screened | 0 / 0 / 0 | 0 |
| 10x Kamil T2 | screened | 0 / 0 / 0 | 0 |
| PacBio T1 | screened | 0 / 0 / 0 | 0 |

No resolved ORF/coding window in these acquired inputs and search bounds. For splice nominations, alternate exon use, intron retention and exonization have not been comprehensively reconstructed. Protein sequence and altered transcript structure remain unknown; reference anatomy above is the available description.
