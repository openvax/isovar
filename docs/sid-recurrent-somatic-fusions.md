# Recurrent somatic-event fusion shortlist

For neoORF follow-up I prioritize the TPST1–CRCP region, FOXO3–STRADA/CCDC47,
and PARD3B–CDKN2B-AS1/CDKN2B. GABBR1–SLC29A1 remains an event watchlist item.
These are likely tumor-associated rearrangements with recurring RNA evidence;
their candidate proteins have not been demonstrated to be translated.

The [published DNA/RNA probe](https://osteosarc.com/fusions/) reports tumor-DNA
support at T1, T1 organoid and T2 for all four event families, with no DNA split
reads in matched blood; PARD3B also has T0 support. It reports RNA junction
evidence at T1, T2 and T3. T1 organoid is a separate library from the same
timepoint, not a third clinical timepoint. Exact ORF comparisons below use
preserved T1 PacBio and T2 tagged ONT inputs; a T3 junction does not establish
the same ORF at T3. Reprocessed alignments are not independent specimens.

| Priority / event | DNA and transcript interpretation | Candidate ORF / altered continuation | T1 PacBio full-ORF fragments | T2 ONT full-ORF fragments |
| --- | --- | --- | ---: | ---: |
| **1. TPST1–CRCP region** | Strong region-level tumor-DNA evidence, but the selected RNA splice has no original DNA VCF call linked to it. TPST1 exon 1 / 5′ UTR joins the CRCP region before the normal TPST1 coding start. | **30 aa + stop; 13-aa shared prefix + 17-aa altered tail**, compared with a potential TPST1 reference uORF. | **15** | **172** |
| **2. FOXO3–STRADA/CCDC47** | Original chr6↔chr17 BND calls in T1, organoid and T2. The representative ORFs start in CCDC47 3′ UTR sequence and continue into FOXO3-region intronic sequence on the opposite strand to FOXO3's coding transcript. | **34 / 36 / 62 / 67 aa + stop**, each with the same **18-aa altered tail** relative to a potential CCDC47 3′-UTR ORF. Other frames remain possible. | **3** for each | **5 / 5 / 3 / 3**, respectively |
| **3. PARD3B–CDKN2B-AS1/CDKN2B** | Strongest longitudinal DNA event: chr2↔chr9 junction in T0, T1, organoid and T2. Intronic PARD3B sequence joins the CDKN2B-AS1/CDKN2B region. The candidate ATG is in the junction insertion. | **10 aa + stop**; insertion-start ORF `MNIRNVHFIS*`. No annotated wild-type protein counterpart is established. | **5** | **11** |
| **4. GABBR1–SLC29A1** | Original chr6 INV-labelled adjacency in T1, organoid and T2. GABBR1-side start is intronic; the nominal partner endpoint has no overlapping supplied Ensembl 115 transcript. | **37 aa + stop** directly linked to the breakpoint. Six codons lie wholly before the join; 31 residues touch or follow it. A reference protein comparison and exact ORF recurrence remain unresolved. | **Unresolved** | **3** |

For the two reference-comparable ORF families, the altered continuations are:

```text
TPST1-region: GRMPILRFSVTTRYFSY*   (17 aa)
FOXO3-region: RYISITEYPYNNGILLSN*  (18 aa)
```

“Altered continuation” measures the suffix after the first difference from a
potential reference ORF. It does not mean every residue in that suffix is an
individual substitution, or that the reference ORF is normally translated.
The TPST1 start is an unannotated ATG 141 nt before the annotated coding start.
The FOXO3-region hypotheses are downstream-UTR ORFs, not an established
FOXO3–STRADA coding fusion protein. Start usage, mature transcript inclusion,
RNA polarity and translation remain unresolved. FOXO3-region and PARD3B ORFs
are supported in the reverse complement of original query sequences.

The 17-aa TPST1 continuation creates 17 eight-mer windows absent from the
supplied potential reference ORFs. The 18-aa FOXO3-region continuation creates
18 such windows, shared by its nested starts. All of those windows are also
absent from the scanned Ensembl 115 annotated proteome. This is a sequence
comparison, not an HLA prediction or proof of tumor-specific antigen display;
normal unannotated translation is not excluded.

## Alternative frames and excluded local ORFs

The FOXO3-region event has seven stopped, breakpoint-linked hypotheses with
at least two full-interval witnesses in each of the two libraries. The table
above emphasizes its nested 18-aa-tail family; alternative frame results are:

| Alternative ORF | Altered continuation | T1 PacBio fragments | T2 ONT fragments |
| --- | --- | ---: | ---: |
| `MAVCRCFSGIFL*` — 12 aa | `GIFL*` — 4-aa divergent continuation of a potential CCDC47 3′-UTR ORF | 3 | 5 |
| `MFFPSWLLQPGGCSWRSSVLYIMAVCRCFSGIFL*` — 34 aa | Same `GIFL*` continuation | 3 | 4 |
| `MPLLYGYSVIEIYRRSNGTQP*` — 21 aa | Unknown; no comparable supplied reference ORF | 3 | 5 |

These frames and nested starts reuse witnesses; their counts must not be
added. Geometry aliases for PARD3B also reuse witnesses and were collapsed.
The 67-aa GABBR1-region hypothesis with nine T1 PacBio fragments is a local
SLC29A1-region path, not a direct bridge across GABBR1's nominal breakpoint.
It therefore does not establish a recurrent GABBR1–SLC29A1 protein.

FTH1–FTL remains an interesting separate 179-aa coding-fusion RNA lead, with
two complete T2 witnesses. It currently lacks patient-level T1/T3 recurrence
of that ORF and linked somatic DNA evidence. IGKC/TRAC signals retain plausible
normal immune-cell origins. Recurrent TALAM1 and repeat/paralog partners need
stronger origin and placement evidence. DSCAML1–ESRRG and SPRED1 have recurring
DNA evidence without positive recurrent RNA evidence in the published probe.

The support counts are complete-interval fragment witnesses within recovered
inputs, not TPM, verified independent RNA molecules or percentages of tumor
cells. PacBio base qualities are unavailable, and recovered inputs may omit
source records. The strict recurring-ORF selection is a positive-evidence
shortlist; a missing entry is not evidence of absent expression.

[The complete curation ledger](data/sid-sv-recovery/curated-somatic-fusions.json)
contains all ten retained exact protein sequences, nucleotide sequences,
reference comparisons, genomic placements, per-library counts, reused evidence
identifiers, selection rules and producer/input hashes. The producer checked
each selected result against the preserved event report's hash. This
supplemental analysis does not change the frozen recovery engine or outcomes.
