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

The shared orfdata catalogue maps positions and REF-validates context SNVs.
Varcode 10.12 adds full-transcript upstream-frame variant consequences and
mapped RNA-product annotation. Use `variant.orf_effects()` for additive DNA
annotations, and the RNA overlay below for reconstructed Isovar candidates.
[Varcode #583](https://github.com/openvax/varcode/issues/583) tracks further
composition and de novo start discovery.

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

## Annotating reconstructed Sid candidates

```python
from isovar import annotate_sv_uorf_evidence, export_sv_rna_orfs
from pyensembl import cached_release

# result is the output of reconstruct_sv_rna, with original mapped RNA paths.
transcript = cached_release(95).transcript_by_id("ENST00000304842")
annotated = annotate_sv_uorf_evidence(result, [transcript])
export = export_sv_rna_orfs(annotated)
```

The overlay validates the entire reference model against the supplied transcript.
Isovar 1.51 uses ORFdata's resolved species/strain, assembly and annotation
context throughout binding and export. The same API works for other species
when supplied a matching catalogue and transcripts. Custom PyEnsembl Genomes
without complete metadata accept an explicit `reference_context=ReferenceContext(...)`
on `annotate_sv_orf_evidence` or the legacy `annotate_sv_uorf_evidence` wrapper.
Conflicting metadata and ambiguous supplied transcript IDs are rejected;
versioned source transcript IDs are looked up and then validated in full.
It converts Isovar's zero-based mapped blocks to Varcode's one-based ancestry,
keeps every candidate's full-interval witnesses, and attaches
`reference_uorf_annotations` to each candidate/occurrence. The export retains
`uorf_catalogue_metadata`, and TSV includes `reference_uorf_annotations_json`.
Absent models and unavailable transcripts are explicit. Empty reference lookup
is not evidence against a sequence-only candidate.

Reproduce the original-read example from the repository with
`python -m examples.osteosarc_sv_orfs output.json --uorf-evidence` (installed
Ensembl 95). The 30/61/72-aa TPST1–CRCP products have 15/1/1 full-interval
witnesses in this packaged PacBio subset. None retains native `KIILSKISV`.
Their RNA polarity, initiation, mutant translation and presentation remain
unresolved. These counts are specific to the subset, not whole-library abundance.
