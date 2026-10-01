# Optional base-quality evidence (#8)

## Scope and invariants

Preserve reported base qualities through allele extraction, read trimming,
transcript classification and overlapping-mate handling. One shared assessment
policy applies to Illumina, ONT, PacBio and unspecified inputs. It accepts an
explicit minimum reported Phred value; it does not choose platform defaults,
calibrate probabilities or reinterpret MAPQ/read-level accuracy as base quality.

Baseline sequence reconstruction, allele counts and cell attribution remain
unchanged, including when qualities are missing. Optional quality assessments
and qualifying support counts are exported separately. A missing score is
unknown, never an invented zero or a passing score. Known low quality fails a
requested threshold; incomplete quality without a known failure is unassessed.

## Evidence representation

- Retain an optional quality vector aligned to each AlleleRead's sequence.
- Keep immutable, original per-alignment allele-quality footprints before
  trimming and merging: original and canonical query intervals, original base
  indices/values, normalized allele and the source alignment identity.
- Repeat normalization must not substitute the canonical position's quality
  for the actual inserted bases. Empty alleles use observed flanking anchors,
  explicitly labelled as anchor quality rather than deletion confidence.
- Mates retain their own footprints. Report individual-read quality outcomes;
  do not multiply probabilities or manufacture independent molecules.
- Metadata never changes observation identity or default output columns.

## API and output

A shared BaseQualityPolicy assesses preserved original allele evidence. Its
optional minimum is descriptive/opt-in; raw supports remain separate from
passed, failed and unassessed support. Integrate it into allele-count and
protein-hypothesis exports and the per-cell reconstruction export, retaining
scope, original evidence sets and unknown-quality counts. Candidate attribution
quality initially qualifies the focal allele only; flanking per-base qualities
are preserved for later calibrated discrimination and are not silently treated
as a whole-ORF confidence score.

## Verification

Use independent CIGAR/quality ledgers for the existing Illumina, ONT and PacBio
fixtures. Test threshold boundaries, missing/partial/constant qualities, mixed
sources, reference/alternate/other alleles, substitutions, original versus
normalized indels, terminal anchors, both strands, trimmed reads, agreeing and
conflicting mates, metadata cloning, legacy objects and default output parity.
Stripping qualities must preserve baseline sequence and cell evidence in
nonconflicting reads. Existing quality-based mate conflict resolution remains
explicit and unchanged; it cannot be reproduced after removing its evidence.

## Primary sources

[SAM specification](https://samtools.github.io/hts-specs/SAMv1.pdf): QUAL is
per-base Phred quality; missing QUAL and MAPQ are distinct concepts.
[PacBio BAM specification](https://github.com/PacificBiosciences/PacBioFileFormats/blob/13.0/BAM.rst#qual):
QUAL and read-level rq have different scopes.
[Dorado Q-score documentation](https://software-docs.nanoporetech.com/dorado/latest/basecaller/qscore/):
base qualities and reported aggregate read quality are distinct summaries.
This PR introduces neither a universal Q20 policy nor calibrated expression
probabilities; those require suitable controls (#411) and weighted-model work (#26).
