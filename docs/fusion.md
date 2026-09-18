# Supplied fusion RNA

Isovar 1.15 adds an explicit sequence-resolved path, separate from SNV/indel
reconstruction. A gene pair, symbolic VCF allele, or DNA breakpoint is not a
fusion transcript. `isovar fusion` does not discover fusions, pad sequences
from reference, select the longest ORF, or assemble a whole transcript.

See [library responsibilities](library-responsibilities.md) for the boundary
between Varcode's hypotheses, Isovar's RNA evidence, and Vaxrank's candidate
evaluation. Event-directed path reconstruction from a BAM is
[`isovar sv-rna`](sv-rna.md); its exploratory candidates are not this validated
input. Reconciliation with Varcode hypotheses remains
[#305](https://github.com/openvax/isovar/issues/305). This guide covers the
supplied-sequence workflow.

```sh
isovar fusion --input fusion.input.json --output fusion.result.json
```

Add `--plot-dir figures/fusions` to save white-background individual PNG/SVG
panels and a vector PDF in a fresh UTC-stamped directory (`--dpi` defaults to
600). Junction RNA and local donor/acceptor annotation are separate panels;
only justified translations receive protein panels. Full original inputs
and model details should accompany published figures.

The JSON input has three keys: `fusion`, `references`, and `reads`. These map
directly to the public `FusionTranscript`, `FusionReference`, and `FusionRead`
dataclasses. Nested mappings use `FusionBreakpoint` and `FusionBlock`.
`reconstruct_fusion(fusion, references, reads)` returns the same serializable
result as the CLI. Defaults are shared: two distinct directly spanning
fragments, and junction peptide lengths 8, 9, 10, 11.

## Coordinates and evidence

All coordinates are plain integers, **0-based interbase**. Genomic intervals
are ascending and half-open, including on the minus strand. Query/cDNA/CDS
coordinates follow the supplied RNA's 5'-to-3' orientation. `reference_name`
identifies the assembly; `contig` identifies the chromosome, without alias
guessing. A breakpoint is the retained partner's boundary, not a VCF POS.

`junction_start:junction_end` is the inserted RNA interval (empty for a direct
join). Blocks must reach both boundaries with the declared chromosome and
orientation. Reference exons are ascending genomic intervals; reference CDS
end includes the stop codon. Supply complete, versioned reference transcript
sequences and all candidate transcript models, not just the preferred one.
Reference-set completeness is the caller's responsibility.

Reads carry sample/library/read/fragment IDs, source, original query offset,
observed sequence, cDNA placement, alignment blocks, and a `mate_number`.
Use `1`/`2` for paired SAM ends that share a query name and `0` for an
unpaired or unknown end. Supplementary/SA mappings must be resolved into one
consistent observation of one mate by the caller; alternative placements must
not be unioned. Paired observations share a fragment ID; processed copies
retain the same library/read/mate identity.
Validated cell/UMI identities may be used as fragment IDs. Isovar validates
exact sequence and mapping agreement, deduplicates processed copies, rejects
conflicting placements, and counts directly spanning fragments separately.
Every supplied cDNA base must have observed read coverage before translation.
Coverage alone is not proof of long-range phase: the supplied assembler/caller
is responsible for that contig hypothesis, which remains in provenance.

Required provenance fields are `sample_id`, `method`, `version`, `parameters`
(including evidence filters), `source`, and `contig_id`. Additional provenance
is preserved. Results retain the sequence digest, event and mapping, read
observations, support counts, annotation identities, and effective parameters.
Read-source authenticity and alignment quality are not inferred from a JSON
claim; retain auditable original alignments and acquisition hashes.

## Frame and result states

Frame transfer is deliberately conservative: the observed donor must exactly
match a collinear interval of an annotated transcript with a complete CDS.
Upstream indels, unknown starts, novel donor paths, and sequence discrepancies
remain unresolved. Nuclear NCBI genetic code 1 is currently supported; other
codes are rejected, not silently translated with code 1.

All compatible donor and acceptor models are retained. Acceptor models label
in-frame versus frameshift context; they do not replace the donor frame.
Actual supplied RNA is translated through its first stop or observed end.
An observed full start has `complete_5prime=true`. Otherwise `cds_start=null`,
and the partial protein is **conditional on the annotated upstream frame**;
the result explicitly records that assumption. `translation_start` is the
first complete codon in the supplied cDNA. `junction_in_translated_cds` is
relative to that start, not to an unobserved full-length CDS. Stop status and
trailing partial codon bases distinguish complete and truncated 3' context.

- `translated`: one translation hypothesis across the compatible supplied
  donor models. This is not a uniquely selected transcript or proven protein
  expression; check partial-CDS and frame evidence before downstream use.
- `ambiguous`: different translation hypotheses or coding/noncoding donor
  alternatives. Downstream tools must not silently rank these.
- `unresolved_frame`: sufficient local RNA support but no justified coding
  interpretation. This is not evidence that no fusion RNA exists.
- `insufficient_support`: absent direct junction support or unobserved cDNA
  bases. A supplied contig alone does not become a validated protein.

Junction peptides strictly cross at least one junction boundary at nucleotide
resolution, including within-codon joins. Insertions have two boundaries.
`downstream_frameshift_peptides` are separate, wholly downstream intervals
whose frame differs from the named acceptor model. Neither list establishes
absence from the reference proteome, antigen presentation, or immunogenicity.

## Real RNA examples and limitations

The pinned osteosarc examples retain original ONT records and exact repeated
junction windows, not reference-concatenated transcripts. TPST1–CRCP has a
donor 5' UTR junction; FOXO3 and PARD3B local windows are intronic relative to
the supplied Ensembl 87 donor models. These remain `unresolved_frame` rather
than acquiring a fictitious coding fusion. Counts describe the selected
exact sequence group, not total event abundance. The PARD3B window preserves
the observed 12-nt inserted sequence. These examples do not establish a
long-read-only detection advantage.

This is additive; existing small-variant APIs/CLI defaults are unchanged.
Fusion results must not be coerced into a fake single-locus `Variant` for
Vaxrank. Symbolic SVs are rejected by the small-variant path (1.14.1).

Sources: [SAM alignment format](https://samtools.github.io/hts-specs/SAMv1.pdf),
[NCBI translation tables](https://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi),
[osteosarc fusion evidence](https://osteosarc.com/fusions/), and
[CTAT-LR-fusion](https://github.com/TrinityCTAT/CTAT-LR-fusion/wiki).
