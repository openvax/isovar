# Matched-normal evidence scope (#595)

Build one typed normal-evidence API shared by the cohort command and acquisition
receipts. Keep an exhaustive source inventory separate from the selected evidence.
Source IDs, pinned snapshots, explicit VCF sample columns, specimen/library IDs,
caller, processing relationships, selection and exclusion reasons define scope;
filename spelling is not a selection rule. Missing assignments stay unresolved.

For each allele, retain per-source records and a reconciled state: observed ALT,
conflicting calls, confidently reference, uncalled, variant-only absence or
unavailable. Resolve source contigs with Varcode once per header. Reference
confidence requires a called reference genotype, accepted record/sample filters,
declared minimum genotype quality/depth and coverage of the whole tested interval.
gVCF END blocks, no-call genotypes, partial calls, other ALT alleles and missing
contigs must not be collapsed into reference. Caller products are not replicates.

Only genotype-indexed literal alleles can enter the supplied background. Keep
caller/filter/genotype provenance and all conflicts; never add every ALT at a
multiallelic site just because one ALT was called. Cohort receipts incorporate
normal inventory/policy digests, nomination evidence and background provenance;
changed scope requires a new analysis directory, never mutation of old receipts.

Audit all VCF products in the pinned October 27 Osteosarc inventory. Inspect
published sample relationships and VCF headers without downloading DNA reads.
Acquire suitable normal call/reference-confidence products or explicitly record
unavailable products and unverified assignments. Generate a versioned Sid normal
assessment and compare its background/origin flags with the original selected union.

Reconciliation exposed a scope regression in nomination-only flanks: distant
spliced exons would lose previously supplied background alleles. Add explicit
background regions to the shared API, pinned manifest and acquisition workflow.
For Sid, use the union of Ensembl 87/115 exons, pinned ncORF coding blocks, previous
background allele spans and nomination flanks. Collect background records once per
source, independently of focal query outcomes; only genotype-indexed accepted
literal alleles are eligible. Region overlap selects complete observed records;
it does not assert reference confidence across their entire REF span. Preserve
the completed v2 receipt and produce a new v3 analysis with expanded scope.

References: VCF specification GT allele indexes and missing calls
(https://samtools.github.io/hts-specs/VCFv4.3.pdf); GATK reference-confidence blocks
and GQ bands (https://gatk.broadinstitute.org/hc/en-us/articles/360035531812-GVCF-Genomic-Variant-Call-Format).

Tests: source labels independent of filenames, stable IDs and processing copies,
MT aliases, alt/decoy exact names, gVCF/no-call/partial blocks, ambiguous contigs,
multiallelic genotypes, disagreement, missing confidence fields, CLI/API equality,
receipt invalidation and original pinned normal-source cases.
