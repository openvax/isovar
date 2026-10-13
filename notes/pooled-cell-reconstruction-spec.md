# Pooled reconstruction and cell attribution (#226)

## Contract

Discover RNA-derived protein candidates with evidence pooled at a declared
sample, library, cell or externally supplied cell-group scope. Keep existing
sample-level results unchanged. An opt-in reconstruction export records each
pool's settings, candidates and no-result reason, followed by attribution of
original observations to the union of candidates. Preserve alternative RNA
contexts, library-scoped cell identities and read provenance.

Discovery's per-base coverage floor is independent of attribution: one read
can support a candidate in a cell even when that cell alone cannot reconstruct
it. Compare every observation against alternative candidates and a reference
allele competitor. Report exact/noisy compatibility, ambiguity, competing
reference evidence, insufficient coverage and unresolved labels explicitly.
Scores are edit costs and compatibility decisions, not calibrated probabilities
of expression. Coverage of a local protein window does not prove a full ORF or
translation in the cell. Counts remain reads, fragments and supplied UMI labels.

## First implementation

- Reuse the existing transcript-aware assembler and its absolute coverage floor;
  preserve its sequence branches rather than taking a majority across cells.
- Add explicit pooling modes and externally supplied, scoped cell groups.
  Unresolved barcodes are accounted for; they cannot create a labelled cell.
- Attribute against the union of discovered nucleotide translations, allowing
  bounded sequence differences and retaining ambiguity between alternatives.
  Report the actually compared interval separately from complete coverage.
- Add an opt-in Python/CLI export with stable candidate, pool and cell IDs.
- Validate pooled single-read cells, group separation, shared sequence,
  noisy reads, reference competition, library collisions, missing labels,
  repeated UMIs, both strands, empty results and original single-cell ONT data.

Automatic biological cell clustering, calibrated cell-expression posteriors,
UMI correction, a new noisy-read consensus assembler, SV reconstruction and
full-transcript discovery remain outside this first PR. These are separate
extensions to #226, not implied capabilities of edit-distance attribution.

## Scientific basis

Pooled transcript discovery followed by per-cell quantification follows the
separation in [IsoQuant](https://ablab.github.io/IsoQuant/cmd.html) and
[FLAMES](https://doi.org/10.1186/s13059-021-02525-6).
[Longcell](https://doi.org/10.1038/s41467-025-60902-2) demonstrates why UMI error,
read truncation and intra-cell isoform diversity require explicit treatment;
matching supplied UMI labels alone does not establish independent molecules.
