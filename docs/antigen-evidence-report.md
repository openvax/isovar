# Antigen-source evidence report

`isovar antigen-report` combines existing Isovar exports for one sample into a
PDF and a machine-readable JSON inventory. It keeps literal DNA alleles,
nominated DNA structural variants, supplied fusion transcripts, selected viral
contigs and targeted antigen genes in separate categories. Their counts have
different denominators and are never summed as “total antigens.”

Install `isovar[report]`. First run the source-specific commands with the same
`--sample-id` and retain their versioned JSON outputs. For example:

```sh
isovar evidence-report --vcf somatic.vcf.gz --bam tumor-rna.bam \
    --sample-id tumor-1 --output small-variants.pdf
isovar antigen-report --sample-id tumor-1 --output tumor-1-antigens.pdf \
    --small-variant-tsv small-variants.tsv \
    --sv-result nominated-sv.json \
    --fusion-result supplied-fusion.json \
    --viral-result viral-rna.json \
    --isoforms-result targeted-antigen-isoforms.json
```

Repeat any JSON input option for additional events or panels. Every supplied
JSON must have the expected versioned Isovar schema and the same `sample_id`.
The command writes `tumor-1-antigens.pdf` and `tumor-1-antigens.json` by
default; `--summary-json PATH` changes the latter path. Omitting an input marks
that category **not assessed**. A zero within a supplied export means the
selected source was assessed and had zero support by that export's criterion.

The small-variant page groups exon-involved alleles by SNV, insertion,
deletion and replacement. It shows original RNA ALT **reads** per allele,
including zero-support alleles and then only positive-support alleles. The
underlying `isovar evidence-report` command also writes these figures directly
from a VCF and indexed RNA BAM/CRAM. A changed reference base meeting any
annotated exon counts; insertions use both flanks as candidate exon involvement.

SV and fusion sections report event-linked RNA junction **fragments** and
separately count translated paths or candidates. They do not reinterpret a
breakpoint-screen split alignment as a small-variant ALT read. The
[HCC1395/COLO829 source PDFs](../reports/hcc1395_colo829/PDF_REPORTS.md) add
an exon-involved SV figure from pinned whole-run breakpoint screens: full
intrachromosomal span for DEL/DUP/INV and nominated breakpoints for INS/TRA/BND.
Those panels show split-alignment candidate counts, including and excluding
zero; they are not direct junction support. The generic composition command
does not compute exon overlap for SV exports because `isovar sv-rna` JSON
does not contain a genomic exon annotation for a complete DNA SV call set.

Viral plots are on separate pages and count selected viral contigs by aligned
and sequence-assigned original RNA **fragments**, each with and without zero
support. They do not use a human exon filter. Selected contigs are the panel
denominator; a panel is not a screen of every possible virus. Coding
predictions are separate from observed protein expression. Targeted
antigen-gene isoforms are reported by nominated gene and RNA-compatible path;
reference-only paths are not counted as observed tumor products. Reference
ncORF overlaps are candidates, not mutant ncORF translations.

The JSON preserves per-source summaries and paths. For large call sets, consult
the original per-call TSV and source-specific JSON for individual events. The
report does not infer somatic status or harmonize cultures, assays, and
reference assemblies; select compatible inputs before composing it.
