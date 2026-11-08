# Shared event and protein semantics (#665, #666)

## Contract

SV joins must follow the nominated adjacency's direction. Ordinary forward
splicing inside a tandem duplication is regional context, not evidence crossing
its backward join. CIGAR deletions may identify a deletion adjacency at its
breakpoints, but incidental alignment indels must neither seed event paths nor
label their proteins as event-linked. Reconstruction, assessment, comparisons,
exports and reports use the same named relation sets. Splice-ambiguous hypotheses
remain exploratory and are not counted as event-linked. Audit summaries retain
exact-breakpoint counts and explicitly report spliced event-compatible counts.

In-frame protein comparisons place the single allele-length-derived gap within
the declared focal interval to minimize substitutions; equal placements prefer
the C-terminal position. Flanking sequence and frameshift tails remain anchored.
Missing/ambiguous residues cannot justify a placement or imply sequence identity.
All consumers use these shared columns/comparisons, including plot highlighting.
The focal nucleotide/codon interval remains provenance, distinct from residues
actually changed. This is a bounded comparison, not unrestricted homology
alignment or a new HGVS annotator.

DNA effect categories are computed from the nominated allele and Varcode effect,
independently of RNA coverage or reconstructed proteins. Reuse one annotation
summary in cohort and evidence-report; keep RNA support as its existing axis.

## Implementation and cleanup

Preserve and incorporate the two existing #665 WIP commits in an isolated branch
based on current master. Consolidate duplicate relation checks and effect rules
in their existing owning modules. Remove obsolete aliases, fallbacks, or derived
sidecar logic superseded by the shared implementation; do not delete unrelated
research artifacts or provenance receipts. Update interpretation documentation
and bump the version. Old scientific checkpoints must not be silently relabelled.

## Evidence and validation

Primary sources checked before implementation:

- [Newman et al., 2015, duplication breakpoint structures](https://pmc.ncbi.nlm.nih.gov/articles/4320257/):
  direct tandem copies retain ordinary within-copy sequence; breakpoint
  orientation determines the novel adjacency. Rejecting a forward join as
  evidence of a backward adjacency is the geometric inference used here.
- [HGVS protein deletions](https://hgvs-nomenclature.org/stable/recommendations/protein/deletion/):
  compare protein sequences and prefer the most C-terminal equivalent position.
  Apply this only inside the existing focal bounds, preserving upstream/downstream
  differences and missing context.

Regressions cover both strands, genuine duplication/deletion junctions, ordinary
splices and 1–3 bp indels, ambiguous vs linked protein counts, H1-2 deletion,
insertion/deletion boundary codons, repeats, unknown context, frameshifts, and DNA
categories with/without ALT RNA. Re-run available Sid examples. Run lint, the
full tests and GitHub CI; merge and deploy from clean master; verify PyPI.

## Replay observations

The original review BAMs retain 96 ONT and 20 PacBio fragments supporting the
TPST1–CRCP adjacency. Removing split partners leaves no event-crossing junctions,
no event-linked annotated proteins and no exploratory event ORFs. The three
reported incidental-deletion ORFs are absent. The H1-2 ONT window reports exactly
five deleted residues (offsets 18–22), no substitutions, and a neutral boundary V.

Streamed support previously filtered indels separately from the shared path
builder. That duplicate filter is removed; final junction filtering now uses
the same event relation policy for both bounded and streamed reconstruction.
No additional evidence sidecars are introduced.

Full-suite validation also guards sequence-assembled joins without a directly
crossing read. Their junction kind is unknown, not a CIGAR deletion; they retain
event linkage with zero direct support. Only explicit D-only evidence invokes
the deletion-placement rule.

The pinned HCC1395 chr9 DUP report disables regional seeding. Its former third
breakpoint-containing path was reached only from an ordinary forward-splice
seed labelled splice-ambiguous. The corrected policy retains the two paths
seeded at the actual adjacency and all seven direct junction fragments. Refresh
the generated junction audit and transcript report to remove the old path count,
FASTA entry and path-limit claim. These are bounded search results, not a claim
that the third sequence is biologically impossible.
