# Optional reported base-quality evidence

Isovar preserves reported per-base qualities from Illumina, ONT, PacBio and
unspecified inputs through RNA allele collection. Reconstruction uses the same
pipeline whether qualities are present or missing. There is no default base
quality threshold or platform-specific cutoff.

`BaseQualityPolicy` adds a separate assessment of the **observed focal allele**.
It does not filter reconstruction, alter raw support or rank protein hypotheses.
It does not score the whole ORF or estimate a cell's expression probability.

```python
from isovar import BaseQualityPolicy, export_protein_hypotheses

export = export_protein_hypotheses(
    results, sample_id="sample-T1", source="rna.bam",
    base_quality_policy=BaseQualityPolicy(min_base_quality=20),
)
```

The threshold here is an example, not a recommended cross-platform value.
`reconstruct_cell_groups` and `CellUmiAlleles.evidence` accept the same optional
policy. Omit it to retain the existing output structure.

The CLI offers the same assessment on allele counts and protein/cell exports:

```sh
isovar allele-counts --vcf variants.vcf --bam rna.bam \
    --min-base-quality 20 --output counts.csv

isovar protein-hypotheses --vcf variants.vcf --bam rna.bam \
    --sample-id sample-T1 --cell-reconstruction sample \
    --min-base-quality 20 --output hypotheses.json
```

## What the categories mean

| Category | Meaning |
| --- | --- |
| `passed` | Every assessed original base has a reported quality at least the requested threshold. |
| `failed` | At least one assessed original base has a reported quality below the threshold. |
| `unassessed` | Quality, an anchor or provenance is missing, without a known failure; or the source mate's allele differs from the merged allele. |

Q0 is a reported value. Missing QUAL, missing entries and BAM's 255 sentinel
are unknown. A partial observation containing both a known failure and a
missing score fails the threshold. Known scores are not filled in from MAPQ,
PacBio `rq`, other reads, or a platform name. Constant reported qualities are
preserved as reported; passing does not establish that they are calibrated.

For an insertion, the assessment uses the **actually inserted bases** in the
original CIGAR, including after repeat normalization shifts the logical allele
interval. For an empty allele (a deletion or reference support at an insertion),
only mapped flanking anchor qualities can be assessed. This is explicitly
labelled `mapped_flanking_anchor_quality_only`, not deletion confidence. An
absent anchor remains unknown; a nearby soft clip cannot replace it.

Overlapping mates keep separate footprints. The existing sequence consensus
still prefers the higher-quality base when mates disagree. An original mate
with a different focal allele is unassessed for the consensus allele; it never
borrows the winning mate's score. Alternative placements of a segment count once:
any known quality failure fails; otherwise every placement must pass to qualify.

## Outputs and counting

Protein-hypothesis and cell-reconstruction JSON add `base_quality_policy` at
the top level and `base_quality` beside raw RNA support. Its `passed`, `failed`
and `unassessed` entries contain RNA support records. Where raw support has an
evidence set, category evidence sets retain hashed, sample/source-scoped read
IDs. Legacy objects without original alignment provenance remain unassessed;
they never receive invented read identities.

Cell candidates and protein attribution also report per-category unique and
ambiguous cell IDs/counts under `base_quality`. These qualify original focal
allele observations already compatible with that candidate. A one-read cell
can appear even when pooled reconstruction needed several reads. Its observed
quality cannot establish the correctness of unobserved ORF bases or distinguish
candidates using only unassessed flanking differences.

Read categories partition the raw support within each record. A fragment, UMI
or cell can have reads in several categories, so those counts are **not
additive**. Raw unique/ambiguous attribution remains visible for cells with
missing quality. Raw protein-hypothesis TSV columns are unchanged; assessments
are in JSON.

Allele-count CSV adds `min_base_quality` and columns such as
`alt_base_quality_passed_reads`, `alt_base_quality_failed_fragments` and
`alt_base_quality_unassessed_reads`. With `--cell-umi-labels --sample-id ...`,
each category also has cell/UMI counts and label completeness fields. Without
the quality option, existing columns and streaming behavior are unchanged.

## Preserved evidence and validation

`AlleleRead.quality_scores` follows its retained sequence. Optional
`source_allele_qualities` contains immutable `AlleleQualityFootprint` records
with original SAM alignment identity, normalized allele, original and canonical
query intervals, original quality indices/scores and base-versus-anchor role.
Those indices always refer to original SAM SEQ/QUAL orientation, including
reverse alignments; they are not offsets into a trimmed or merged read.
`BaseQualityPolicy.assess(read)` exposes the per-source assessments.
Quality metadata does not change observation equality or default read tables.

Offline tests use independently recorded CIGAR/quality ledgers from public
Illumina/STAR, ONT/minimap2 and PacBio/minimap2 reads, including a real insertion
whose original Q35 would incorrectly become Q28 if canonical coordinates were
used. Tests also remove QUAL, compare unchanged alleles and reconstructed ORFs,
and verify sparse cell attribution, both strands, trimming, indels and mates.
Removing QUAL cannot preserve an existing quality-dependent mate conflict
resolution when the deciding evidence itself is removed.

The [SAM specification](https://samtools.github.io/hts-specs/SAMv1.pdf),
[PacBio BAM specification](https://github.com/PacificBiosciences/PacBioFileFormats/blob/13.0/BAM.rst#qual)
and [Dorado Q-score documentation](https://software-docs.nanoporetech.com/dorado/latest/basecaller/qscore/)
distinguish reported per-base qualities from alignment or aggregate read
quality. Automatic calibration requires validated controls
([#411](https://github.com/openvax/isovar/issues/411)); quality-weighted candidate
likelihoods and whole-ORF confidence remain future work.
