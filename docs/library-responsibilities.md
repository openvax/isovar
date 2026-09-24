# Varcode, Isovar, and Vaxrank

Three OpenVax libraries share the work of turning tumor sequencing into vaccine
candidates. Isovar owns the RNA evidence: it reconstructs the sequences the
reads actually show and works out what the reads can and cannot distinguish.
An effect label from annotation alone cannot answer those questions.

| Library | Responsibility |
|---|---|
| **Varcode** | Generate structural/transcript hypotheses and predict their coding consequences. |
| **Isovar** | Reconstruct RNA-supported sequences, compare them with those hypotheses, and preserve unresolved alternatives. |
| **Vaxrank** | Evaluate the resulting protein/peptide candidates, retaining their evidence. |

## How RNA reconstruction should work

These steps describe the complete workflow for SVs, which is tracked in
[#305](https://github.com/openvax/isovar/issues/305). `run_isovar` does not
provide it; [`isovar sv-rna`](sv-rna.md) is its first part, and
[Available today](#available-today) lists what exists now.

1. **Collect anchored evidence.** Include relevant aligned sequence, usable soft
   clips, observed supplementary alignments, mates, and long reads. A clipped
   tail alone is not a partner assignment.
2. **Assemble alternative paths.** Preserve supported isoform/haplotype branches;
   do not fill unobserved gaps with reference bases and call them observed RNA.
   Accept external assemblies through the same downstream checks.
3. **Compare competing explanations.** Reconcile assemblies with Varcode's
   hypotheses, normal transcripts, and other genomic matches. Allow a supported
   assembly to introduce a junction absent from the original candidate list.
4. **Check the original reads.** Validate sequence and junction support. Count
   read segments, fragments, and independently established molecules separately;
   supplementary or alternative records are not additional molecules.
5. **Return sequence and uncertainty.** Predict protein where the frame is
   justified. Keep unresolved frame/sequence results and multiple compatible
   proteins, rather than choosing an isoform because it appeared first.

Every result should keep:

- the nucleotide-level structure, reference and annotation, and transcript IDs;
- sample and library identity, and read and fragment provenance;
- support, conflicts, how complete the sequence is, and what the frame assumes.

Grouping identical proteins must keep the structures and evidence behind each.

Four outcomes must stay distinct: *supported*, *contradicted* within an observed
region, *indistinguishable* over the covered region, and *not covered*. A long
read may resolve a branch that shorter reads cannot, and several expressed paths
can coexist. When a resource limit stops a search, the unexamined hypotheses
must be reported as unexamined, not as unsupported.

Technical sequence must not become fusion evidence merely because it aligns
somewhere. Explicit adapter/poly-A annotation and trimming are available
([read-end inference](read-end-inference.md)); recognizing adapter-derived
sequence in reused soft clips remains [#302](https://github.com/openvax/isovar/issues/302),
and inferring an unknown kit remains [#309](https://github.com/openvax/isovar/issues/309). High base quality does not
establish biological origin, and missing qualities do not by themselves reject
an otherwise usable observation.

## Available today

- `run_isovar` reconstructs local RNA context for small variants and keeps
  protein alternatives in `IsovarResult.sorted_protein_sequences`.
  `top_protein_sequence` is a selection convenience, not a unique biological answer.
  The small-variant read path rejects symbolic SVs.
- [Nominated-SV RNA paths](sv-rna.md) use `reconstruct_sv_rna` / `isovar sv-rna`.
  From a BAM and annotated models, they recover junction-seeded paths,
  including spliced joins whose DNA breakpoint is intronic. They transfer
  annotated frames into novel downstream sequence and label sequence, frame
  and event linkage separately, and `--predictions` compares supplied protein
  predictions with those paths. Full reconciliation with Varcode hypotheses,
  competing genomic placements and Vaxrank consumption remain open.
- [Supplied fusion RNA](fusion.md) uses `reconstruct_fusion` / `isovar fusion`.
  It validates externally supplied sequence, mappings, and read support and
  retains compatible reference/frame interpretations. It does **not** discover
  or assemble a fusion from BAM soft clips. Its result uses the same
  [RNA path format](sv-rna.md#the-rna-path-format) as `isovar sv-rna`, and
  Vaxrank's `fusion_antigens_from_isovar` turns its coding hypotheses into
  candidate antigens without choosing an isoform.
- Matched germline/co-somatic attribution of reconstructed edits remains
  [#297](https://github.com/openvax/isovar/issues/297). Sequence reconstruction
  and identifying the origin of each edit are different tasks.

RNA supports an RNA sequence, not direct protein translation. Local junction
support can justify a partial protein without resolving a full mature transcript.
Transfer structures and evidence to Varcode for consequence prediction; downstream
Vaxrank evaluates the protein/peptide candidates without redoing RNA inference.

## Other library guides

- [Varcode: hypotheses and coding consequences](https://openvax.github.io/varcode/library_roles/).
- [Vaxrank: candidate evaluation and evidence](https://openvax.github.io/vaxrank/library-responsibilities/).
- [Fusion assembly and read-back validation methods](https://pmc.ncbi.nlm.nih.gov/articles/PMC6802306/).

