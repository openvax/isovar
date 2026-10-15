# Recover observed viral RNA from combined-reference alignments

Close the immediate #465 gap: a BAM aligned to GDC GRCh38.d1.vd1 should yield
observed viral nucleotide sequences without a nominated mutation or a human
GTF. Add a Python entry point and `isovar viral` command using the existing
coordinate-aware overlap assembly, disk-backed input index and original-read
support machinery. Do not model an entire viral genome as an RNA transcript.

Package the official GDC viral-contig panel with acquisition checksums and
reconcile its names against the existing published sequence dictionary. Validate
the combined BAM/FASTA dictionary and the selected viral DNA checksums. Default
to the complete panel; allow an explicit contig subset. Custom references must
declare their viral contigs explicitly. Human-only alignments are an error for
the GDC profile, not evidence of viral absence.

Discover exact observed paths with bounded per-contig budgets so abundant reads
on one virus do not prevent another virus from being examined. Count support
from the full declared input after discovery, retaining original read identities,
mates, secondary/supplementary placements and missing-placement ambiguity.
Preserve splice paths and observed alleles; never fill gaps with reference DNA.
Keep counts by contig distinct from counts per reconstructed path, and report
unassigned evidence, filtered records, truncation and zero eligible support.

Export reconstructed nucleotide FASTA plus versioned JSON with reference/source
identity, mapping, coverage, witnesses and completeness limitations. RNA strand,
full transcript ends, viral presence in tumor cells and protein coding are not
established by alignment alone. No implicit six-frame proteins or unannotated
viral gene names. Sequence-matched viral coding annotations and general coding
exceptions remain separate work in #465; ordinary annotated targets can already
use the isoform API.

Validate on offline cases for overlapping reads, splice/allele differences,
opposite placements, human/viral multimapping, paired fragments, excluded reads,
empty results, missing references, checksum mismatches, late support and discovery
limits. Include a pinned real viral reference sequence with independently known
synthetic read truth, and exercise the installed CLI without a GTF. Run lint,
full tests and CI; bump the version, merge and use the standard PyPI deployment.

Primary sources consulted before implementation: GDC reference files and viral
sequence list, and NCBI feature-table documentation. GDC's viral list includes
historical/revised reference accessions; matching names alone cannot establish
sequence identity or supply transcript/CDS annotations.

Validation exposed two shared defects which must be fixed in this PR: known
missing SA/mate records could be treated as unique (#509), and target-oriented
observation indexes did not provide the mirrored views required by 5-prime
extension (#510). Preserve original identities while correcting both; the
actual-GDC experiment must recover every independently specified sequence.
