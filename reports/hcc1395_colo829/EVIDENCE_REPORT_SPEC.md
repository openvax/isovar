# Variant and RNA evidence reports

## Inputs and semantics

`isovar evidence-report --vcf calls.vcf[.gz] --bam tumor.bam --sample-id NAME
--output report.pdf` uses Isovar's read collection and local protein context
translation on a variant set and an
indexed RNA alignment. The report records the selected set, original RNA allele
support and local protein contexts; a sequence observation is not proof of a
complete isoform or protein expression. A per-variant TSV accompanies the PDF.
ReportLab is an optional report dependency.

The HCC1395 and COLO829 source audits continue to rebuild from their committed
whole-run tables. Their annotation table gains an explicit Ensembl 87 exon
overlap flag. It is computed against every exon, regardless of coding status.
SNVs and replacements query changed reference bases; deletions query deleted
reference bases, even when neither endpoint is exonic. An insertion has no
changed reference base, so either immediately flanking base counts as a
candidate exon overlap and is labeled as such. Coding effects and orfdata
coding-block overlaps remain separate annotations.

## Figures and denominators

Every report has a nested funnel: all ingested DNA source records, then
small-variant records, variants with any exon involvement, exon-involved
variants with observed RNA ALT, and those with a local mutant protein context.
When no SV source records are supplied, the two initial stages coincide and
are shown once. Report the absolute count and percentage of the initial input
at every stage. SV source records leave the funnel before exon annotation;
their RNA screen is separate. Additional figures show RNA ALT-read support by
variant class, both with zero-support calls and conditional on positive
support, and the number of changed amino-acid positions in the top-ranked
local protein context per RNA-supported variant, both with and without the
single-residue category. The latter is a local window measure: it cannot
infer an unobserved frameshift tail or full-length isoform. RNA statuses are
mutually exclusive; reference coding, ncORF overlap and RNA support are not.

The HCC1395/COLO829 reports count DNA-predicted frameshifts, reference ncORF
overlap candidates and RNA-resolved SV protein products separately. They show
highest-support and non-SNV examples from the computed tables. A "least
self-like" ranking uses the fraction of mutation-overlapping 9-mers absent
from a pinned reference-protein FASTA. It is labeled as a reference-proteome
exact-match screen, not patient-specific self, HLA binding or immunogenicity.
An unavailable reference FASTA leaves this category unevaluated; no substitute
score is silently invented. The reports retain source, coding, ncORF and SV
tables and caveats.

## Verification

Check exon boundaries and an intron-to-exon deletion with synthetic tests;
check funnel nesting and denominators independently of PDF rendering. Run the
new command on a bounded VCF/BAM fixture, parse its TSV and PDF, regenerate
both published PDFs, and visually inspect every page. Then run repository lint,
tests and CI before shipping.
