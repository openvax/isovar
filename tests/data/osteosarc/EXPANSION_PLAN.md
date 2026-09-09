# Issue #218: complete vaccine-cohort RNA audit

## Contract

Pin every vaccine-included variant in the public index, its exact genomic
alleles, and all available RNA alignment products discovered from the BAM
catalogue and data-page inventories. An alignment product is not necessarily
an independent sample. Preserve assembly, processing and donor/timepoint
ambiguities explicitly. Never query GRCh38 positions in a GRCh37 alignment.

Acquire original indexed regional reads with bounded, resumable downloads,
checksums and source/header provenance. Do not download entire BAMs, alter
records, combine technical reprocessings, or infer patient-specific counts
from an unresolved pooled library. Keep small offline fixtures separate from
the full-region matrix used for counts.

Run Isovar's read-to-ranked-protein path for every resolved variant/source,
under documented default and primary-only settings. Report independent
CIGAR observations and quality sensitivities separately. Preserve explicit
null-count source/errors, no overlap, no callable allele, zero alternate
support, and sequence/matching/translation failures. Counts must name their
units; compatible template names are not full-window or UMI molecule counts.

Validate coding outputs against independent edits of pinned, original
annotation/cDNA/protein sequences, checking strand, frame, mutation interval,
deletion junction and stops. Scope transcript expectations explicitly; use
the appropriate nuclear/mitochondrial genetic code. Distinguish additional
RNA differences from incorrect translation, and reconstruction from final
vaccine selection. Keep insertion/complex stress cases separately labelled
when they are not vaccine members.

## Implementation sequence

1. Snapshot metadata and reconcile source/variant inventories. Record source
   and biological-identity gaps before interpreting any counts.
2. Implement bounded regional acquisition and deterministic provenance;
   independently pin reference subsets and variant edits.
3. Generate the full matrix and readable report; investigate unexpected
   failures, file confirmed bugs, and retain original regression witnesses.
4. Add offline integrity, failure-state, CIGAR and protein-oracle tests plus
   representative real-read composed-pipeline regressions. Reproduce outputs.
5. Bump version, pass lint/full suite and PR CI, then follow the repository's
   merge/deploy workflow. Do not claim corpus completion while required
   sources or protein checks remain silently unexamined.

## Findings that change sequencing

The first full GRCh38 pass exposed #223: MT variants cannot find chrM reads.
Preserve its before-fix matrix, repair contig resolution with focused tests,
then regenerate the audit before judging mitochondrial protein results.
Name aliases alone never justify hg19/rCRS mitochondrial coordinate reuse.

The user's mitochondrial follow-up exposed #224: inherited AGA/AGG arginine
assignments disagree with NCBI table 2 and direct 2023 mtRF1 experiments.
Correct and exhaustively test the table before final protein validation.
The MT-ND5 A220T window itself contains TGA/ATA but no AGA/AGG.

NUMT origin and clonality are separate, unresolved questions. Label the
mitochondrial protein expectation as conditional on mitochondrial origin;
do not turn correct translation, MAPQ, or high RNA allele fraction into
proof of organelle origin, DNA heteroplasmy, or cancer-cell fraction.
Retain read-level origin evidence and known-NUMT overlap as diagnostics,
not an automatic exclusion rule. RNA VAF remains an expression-weighted
observation; mitochondrial copy number/heteroplasmy are not diploid CCF.

The restored high-depth ONT mitochondrial evidence exposed #227: adaptive
candidate/support matching can take tens of minutes and several GB. The
full audit now uses an explicit 120-second budget per source/locus/mode,
retaining completed read counts and the interrupted stage. This is an audit
resource limit, not a production threshold, biological negative, or fix for
the performance issue. Do not downsample full-matrix evidence to hide it.

## Primary semantics

- https://samtools.github.io/hts-specs/SAMv1.pdf
- https://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi
- https://osteosarc.com/variants/
- https://osteosarc.com/bams/bams.json
- https://osteosarc.com/data/
- https://registry.opendata.aws/sid-osteosarc/
