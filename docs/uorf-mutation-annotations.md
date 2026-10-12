# Interpreting mutations in upstream ORFs

An interpretable annotation has four independent layers. A favorable sequence
context is a prior about initiation, rather than an experimental observation.
No context threshold removes an ORF or overrides its ribosome/MS evidence.

## 1. Identity and location

Pin the assembly, reference release, transcript version, ORF ID, spliced coding
blocks and protein sequence. Reference evidence belongs to this model, rather
than every transcript of the gene. Use the actual patient's RNA path when an
edit changes splicing, orientation or fusion partners.

Locations can overlap. Number initiation-context positions with the A of ATG
at +1, without a position zero. A variant at +4 is both a context mutation and
a mutation of the second coding codon. A variant upstream of the start can
change initiation context without changing a protein residue.

## 2. Sequence consequence

| Event | Sequence interpretation | What remains unresolved |
|---|---|---|
| Start-context SNV | Compare every observed preferred position before/after the edit; retain missing positions. | Initiation efficiency and downstream CDS regulation are not measured. |
| Start-codon mutation | Report ATG gained/lost or another codon change. | Loss of ATG does not prove no translation; near-cognate or alternative starts remain possible. |
| Coding SNV/MNV | Translate in the particular uORF frame; distinguish synonymous, missense and stop changes. | A canonical-transcript UTR label does not supply the uORF consequence. |
| Indel | Reconstruct the altered spliced sequence and frame, then follow it to a stop or a sequence boundary. | An anchor alone does not specify the full affected protein interval. |
| Stop loss | Extend along the actual transcript path to the next in-frame stop. | Incomplete downstream sequence means an unresolved extension, not a complete protein. |
| Splice/fusion/SV | Reconstruct every supported RNA path, retained native frame and new sequence. | An initiation prior does not establish transcript inclusion, DNA-event linkage or somatic specificity. |

Isovar's catalogue API currently maps coding positions and REF-validates SNVs
in the retained initiation window. It does not yet reconstruct arbitrary uORF
indels, stop extensions or rearranged proteins from catalogue entries. Those
require full reference/patient transcript reconstruction. The proposed Varcode
implementation is tracked in [Varcode #583](https://github.com/openvax/varcode/issues/583).

## 3. Initiation plausibility

The retained context is up to nine transcript bases before the start, the
three-base proposed start codon, and up to six bases after it. Each preferred
Kozak position has an observed base, a preferred base/set and a true/false/null
match. The optional outer GCC is displayed alongside the core positions.
The two commonly emphasized positions are purine at -3 and G at +4.
The key-position counts are descriptive, not a calibrated efficiency score.
A mismatch remains compatible with initiation, and partial sequence is not a
negative result. Non-ATG codons are retained without estimating their efficiency.

These conventions are informed by [Kozak's mutagenesis experiments (1986)](https://doi.org/10.1016/0092-8674(86)90762-2)
and [vertebrate sequence survey (1987)](https://doi.org/10.1093/nar/15.20.8125).

If the first three supplied bases are ATG, there are no observed upstream
bases. Keep the candidate, mark its context partial, and distinguish a supplied
sequence boundary from a verified biological 5-prime end. Only the latter
establishes cap distance. [Short-leader experiments](https://pubmed.ncbi.nlm.nih.gov/1820208/)
found increased bypass of the first AUG as the leader shortened from 32 to
three bases; they do not measure every possible zero-length-leader transcript.
[Special short-leader mechanisms such as TISU](https://pmc.ncbi.nlm.nih.gov/articles/PMC3177215/)
also show why Kozak matching alone cannot decide initiation. An RNA read or
incomplete assembly starting at ATG does not demonstrate a leaderless mRNA.

## 4. Evidence and abundance

Keep reference Ribo-seq, shotgun MS and HLA MS separate, with paper/data pins,
peptide intervals and source/review uncertainty. Overlay native peptide
coverage only on retained native residues. A peptide in a discarded native
tail cannot support the altered fusion tail. Separately record patient RNA
support, phase, translated mutant peptide evidence and HLA presentation.
This reference catalogue contains no patient RNA abundance estimates.

For example, native TPST1 has `GAGGCCAGG[ATG]CCGTCC`: A at -3, C at +4.
The -3 A>C SNV loses a preferred initiation-context base without changing a
coding residue. The +4 C>G SNV gains the usual preference and changes the
second codon from CCG (P) to GCG (A). Neither comparison measures a change in
initiation. The native MS peptide `KIILSKISV` covers residues 31-39, not these
start-region bases and not the TPST1-CRCP fusion's altered tail.

## Varcode improvement path

Add an assembly/transcript-aware ORF-model provider or neutral interchange
format, with no Varcode dependency on Isovar. Keep canonical effects and add
per-ORF consequences with stable ORF identities. Reuse Varcode's existing
transcript-layout, mutant-transcript and hypothesis machinery for coding,
splice, fusion and phased haplotype changes. Export products, alternatives,
source pins and unresolved reasons together. A sequence-only newly created
start must remain distinguishable from a reference model with experimental
support. See [Varcode #583](https://github.com/openvax/varcode/issues/583) for
staging and acceptance cases.
