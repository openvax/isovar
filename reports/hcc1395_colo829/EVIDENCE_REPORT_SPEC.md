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

The cell-line reports have separate nested funnels for small variants and
structural variants. Each stage reports a count and percentage of its own
DNA input class. The small-variant path is exon involvement, RNA ALT read,
and local protein context. The SV path is exon involvement, RNA
alignment near a nominated breakend, same-read split or cross-end mate
candidate, resolved junction sequence, and translated SV product. Breakend
coverage alone is not event support. The generic VCF/BAM command has only
the small-variant funnel because it does not interpret symbolic SVs.
Additional figures show RNA ALT-read support by
variant class, both with zero-support calls and conditional on positive
support, and the number of Isovar-marked mutation-interval amino-acid positions
in the top-ranked local protein context per RNA-supported variant, both with
and without the single-residue category. The marked interval may be unchanged
at the amino-acid level for a synonymous allele; the figure is not a
reference-verified change count. It cannot infer an unobserved frameshift tail
or full-length isoform. RNA statuses are
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

For selected high RNA-support products, show an Ensembl 87 transcript's
reference local protein window beside Isovar's top mutant window. Align at
the mutation boundary, highlight divergent residues, and distinguish a
frameshift's boundary-anchored comparison from a true downstream homology
alignment. Reference contexts are pinned in a small TSV so PDF regeneration
does not require a live or locally cached reference.

## Verification

Check exon boundaries and an intron-to-exon deletion with synthetic tests;
check funnel nesting and denominators independently of PDF rendering. Run the
new command on a bounded VCF/BAM fixture, parse its TSV and PDF, regenerate
both published PDFs, and visually inspect every page. Then run repository lint,
tests and CI before shipping.
