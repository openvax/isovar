# Shared BAM/FASTA reference validation (#648)

Resolve the FASTA spelling once per BAM contig with Varcode's exact-first
resolver. Preserve the BAM name for RNA coordinates and diagnostics; use the
resolved FASTA name for lengths and target-sequence MD5 checks. Record the
mapping in reference provenance. Reject ambiguous/missing names and a mapping
that collapses distinct BAM contigs onto one FASTA entry.

Review found that validating aliases alone left sequence consumers using the
unresolved BAM spelling. Resolve against the FASTA in transcript extraction
and viral coding-reference loading too; keep BAM names in returned RNA models.
The regression tests exercise the CLI and viral reconstruction, not just the
validation helper. Ship this independently of the CTA catalogue as 1.76.2.

Keep the named GDC profile's exact dictionary requirement. For a custom
reference, equal lengths and accepted aliases remain insufficient to establish
sequence identity: target contigs without BAM M5 remain explicitly unverified.
The SAM specification defines SN, LN and M5 separately:
https://samtools.github.io/hts-specs/SAMv1.pdf.

Test canonical, mitochondrial, accession/UCSC, viral and decoy names, exact
precedence, missing/ambiguous names, and length/checksum failures. Then rerun
Sid's failed preparation without changing the original BAM or its identity.
