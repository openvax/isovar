# Extend the focused Sid RNA subsets to non-SNV candidates

The October 5 focused archive covers vaccine claims and four named fusions,
but "interesting variants" also includes the non-SNV examples documented in
`figures/osteosarc/neo-orfs-and-long-reads.md`. Several vaccine-listed frameshifts
are already present; acquiring their reads did not complete a protein audit.

Add an optional expanded selection containing all indel and splice-site
nominations in the pinned small-variant metadata, plus the exact FOXO3–STRADA/
CCDC47 and ATP5MG–KMT2A breakpoints. Keep every distinct allele, including the
two GLIS3 and MAP2 nominations, and report nonliteral/unresolved alleles as
unassessable rather than omitting them or declaring negative RNA support.

ATP5MG–KMT2A is absent from the frozen catalogue. Preserve that catalogue and
create a derived inventory with the checksum-pinned existing original-RNA
fixture as an explicit supplemental nomination. Record its retained-side
geometry from the fixture's interbase breakpoints, without inferring DNA origin
or assigning a new protein. Keep the parent inventory and fixture provenance.

Acquire only new regional inputs from the four previously selected products,
with original qualities, tags, mates and supplementary records. Reuse the same
bounded recovery and 8 GiB free-space guard. Keep the original archive immutable;
publish a separate addition ledger and bundle with exact source/target
denominators, checksums, native BAM recounts and independent candidate checks.
The complete SV catalogue and general cryptic/intron-retention/exonization ORF
discovery remain separate, pending work. Lint, full tests, CI, merge and PyPI
deployment are required for the follow-up PR.

The user subsequently clarified: "Just look for SVs and neoORFs." Focus the
result on actual RNA-derived candidates. Screen the completed SV reconstructions
and the indel/splice inputs for novel sequence, retain sample-specific evidence,
and compare candidate peptides to checksum-pinned Ensembl 115 proteins. Check
full candidate witnesses independently. Report altered upstream and frameshift
ORFs separately from short unframed hypotheses and splice-compatible joins;
none establishes translation, tumor specificity or antigen presentation.
