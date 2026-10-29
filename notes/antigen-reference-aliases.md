# Shared BAM/FASTA reference validation (#648)

Resolve the FASTA spelling once per BAM contig with Varcode's exact-first
resolver. Preserve the BAM name for RNA coordinates and diagnostics; use the
resolved FASTA name for lengths and target-sequence MD5 checks. Record the
mapping in reference provenance. Reject ambiguous/missing names and a mapping
that collapses distinct BAM contigs onto one FASTA entry.

Keep the named GDC profile's exact dictionary requirement. For a custom
reference, equal lengths and accepted aliases remain insufficient to establish
sequence identity: target contigs without BAM M5 remain explicitly unverified.
The SAM specification defines SN, LN and M5 separately:
https://samtools.github.io/hts-specs/SAMv1.pdf.

Test canonical, mitochondrial, accession/UCSC, viral and decoy names, exact
precedence, missing/ambiguous names, and length/checksum failures. Then rerun
Sid's failed preparation without changing the original BAM or its identity.
