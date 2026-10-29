# Frames on observed RNA paths

`isovar.transcript_path_frames` provides the transcript coordinate and CDS-anchor
layer used by SV reconstruction and path assessment. It can also assess an
already observed intronic RNA path without inventing an SV nomination. Read
acquisition, assembly, allele calling and original-witness coverage remain the
responsibility of the caller.

```python
from isovar.transcript_path_frames import TranscriptPathModel, assess_transcript_path

model = TranscriptPathModel(reference)  # an explicitly versioned FusionReference
assessment = assess_transcript_path(
    observed_sequence,
    observed_positions,  # one (contig, zero-based position, strand) or None per base
    model,
    min_anchor_bases=12,
    permitted_departures=observed_junctions,  # query offset -> assessed relation
    allow_retained_introns=True,
    focal_interval=allele_query_interval,
)
for reading in assessment.readings:
    print(reading.amino_acids, reading.focal_coverage, reading.complete_5prime)
```

Inputs describe the RNA in its 5′ to 3′ orientation; minus-strand genomic
positions therefore decrease within each aligned block. The supplied sequence
requires explicit DNA bases (`ACGT`); ambiguous bases raise before translation,
matching SV reconstruction's existing sequence eligibility. Missing placements
are `None`. All intervals are
zero-based and half-open. Query coordinates never become chromosome coordinates
by stripping a prefix or parsing a contig name. A local-reference alignment
requires models and variant coordinates in that same explicit coordinate system.

`TranscriptPathModel` caches immutable transcript positions, an immutable offset
lookup and the supplied CDS protein. Reuse it across paths. Its input is the
existing validated `FusionReference` type: versioned exons and spliced RNA, with
optional CDS bounds, currently restricted to standard genetic code 1. Noncoding
models keep their missing CDS. No longest-ORF substitution is made. Models and
assessments support pickling; `dataclasses.asdict` exports their declared fields
without copying derived indexes.

`TranscriptFrameAssessment` and `PathFrameReading` are frozen typed results.
The assessor transfers a frame downstream from an exact, collinear coding
anchor. Sequence after that anchor can depart from annotation only at an
explicitly permitted departure, or at an observed exon-to-intron continuation
when `allow_retained_introns=True`. This latter interpretation is specific to
each transcript model: another supplied transcript may explain the same path
as an ordinary exon. Assess every credible model and retain those alternatives.

| Result | Meaning |
|---|---|
| `readings` | Coding-anchor hypotheses whose first positional departure is permitted. |
| `collinear_readings` | Coding-frame hypotheses without a positional departure before their observed end/stop; they can still contain substitutions. |
| `other_departures` | Coding-anchor attempts that departed elsewhere and were withheld. |
| `reason="no_sequence_anchor"` | No exact sequence/placement anchor reaches the requested length. |
| `reason="CDS_unavailable"` | The RNA anchors to this model, but it supplies no CDS. |
| `reason="no_coding_anchor"` | Exact matches exist, but no sufficiently long coding anchor exists. |

An ordinary collinear reading is not automatically reference-identical or a
mutant protein. Use the existing [protein comparison](protein-hypotheses.md)
and provenance rules when constructing an export. Different anchors ending at
the same codon share their longest reading, retaining all anchors. Different
frames remain distinct. The generic assessor does not count supporting reads.

`translation_end` includes an observed stop codon; `amino_acids` excludes it.
`complete_5prime` means the supplied CDS start is present in the observed path,
not that initiation was measured. Otherwise the frame is conditional on an
upstream CDS anchor. No missing RNA bases are filled from reference.

The optional focal interval is an already established allele's location in the
path. A deletion is an empty boundary and needs translated RNA on both sides.
`focal_coverage` is `covered`, `partially_covered`, `upstream_of_frame_anchor`,
`after_stop`, or `beyond_observed_translation`; without the interval it is
`not_assessed`. A reading that stops before the allele is not evidence for an
allele-containing protein. Even `covered` alone does not establish a changed
amino acid, end-to-end molecular phase, translation, or mutation-induced splicing.

The SV adapter retains the existing JSON schema, grouping, permitted-junction
rules and default outcomes. This release adds no CLI mode and does not enable
retention discovery automatically. The next original-read consumer is tracked
in [#583](https://github.com/openvax/isovar/issues/583), following the
[local RNA assessment](variant-read-assessment.md). The distinction between an
observed splice path and allele-associated splicing follows the
[MiSplice primary study](https://www.nature.com/articles/s41467-020-19307-6).
