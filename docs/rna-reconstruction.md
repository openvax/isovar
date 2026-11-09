# Shared RNA reconstruction and support

Isovar separates original RNA evidence, candidate discovery, support measurement,
and biological interpretation. Small variants, structural adjacencies, CTA
isoforms and viral paths reuse the public contracts below. An event adapter
supplies its seeds, references and interpretation rules; it does not define a
new read identity or silently reinterpret another workflow's count.

## Public engine

| Module | Responsibility |
|---|---|
| `read_collector` | Shared eligibility and technical-sequence policy |
| `read_identity`, `rna_evidence` | Original segment/fragment identities, scoped counts and evidence-set unions |
| `rna_library` | Validated original-library indexes and caller/context ownership |
| `rna_observation` | Original aligned observations, alternative placements and bounded geometry reuse |
| `rna_path_geometry` | Query-to-genome placements, collinear runs, overlap routing and geometric spans |
| `rna_path_assembly` | Exact overlap extension with consistent original-segment choices |
| `rna_reconstruction` | Bounded path discovery and complete-input support for fixed candidates |
| `transcript_path_frames` | Conditional coding frames on observed paths, without reference-filled sequence |

The shared engine does not import the SV, CTA or viral adapters. SV acquisition,
adjacency classification and result assembly live in `sv_rna_collection`,
`sv_rna_event` and `sv_rna_result`. `sv_rna` orchestrates those operations.
`sv_rna_witnesses` measures junction-linked interval witnesses separately from
ORF enumeration. Bounded and streamed SV support consume these same contracts.

The regular small-variant route remains allele anchored: `ReadCollector`,
`VariantSequenceCreator` and `ProteinSequenceCreator` enforce the focal allele
and coding-context requirements. It shares original identities and support
accounting with the other routes. Extended literal-allele paths additionally use
`allele_path_seeds` and the same path discovery/support engine as CTA and viral
reconstruction. A small allele is never represented as a fictitious SV merely
to invoke shared code.

## Discover first, measure against the original input

```python
from isovar.read_collector import ReadCollector
from isovar.rna_library import original_library_context
from isovar.rna_reconstruction import discover_rna_paths, support_rna_paths

collector = ReadCollector.for_rna_paths()
regions = {"1": [(1000, 2000, "+")]}

# bam is an open pysam.AlignmentFile; scratch is a disposable directory.
with original_library_context(
    bam, collector, scratch,
    library_index=None,
    reference_name="GRCh38",
    read_filter_identity=None,
) as library:
    paths, discovery = discover_rna_paths(
        library, regions, collector,
        max_discovery_reads=10000, max_paths=256, min_overlap=30,
    )
    candidates = [
        dict(id=f"path-{i}", sequence=sequence, placements=placements, coding=[])
        for i, (sequence, placements) in enumerate(paths)
    ]
    support = support_rna_paths(
        library, candidates, regions, collector, bam.header.to_dict(),
        sample_id="sample-1", source="original-library-1", label="path",
    )
```

The example deliberately supplies no coding model: these are observed RNA paths,
not automatically protein antigens. Callers attach sequence-matched references
and coding interpretations through the relevant event adapter.

`discover_rna_paths` returns paths and explicit budgets, considered/selected
segments, rejection counts and truncation. `support_rna_paths` takes those fixed
candidates, visits the complete declared input and updates each candidate with
its measurements. Its returned dictionary contains target/unassigned fragment
totals, ambiguity classes and original witnesses. Discovery truncation does not
truncate the later support pass, but can leave undiscovered alternatives.

The `label` argument preserves native schema spellings such as `isoform_ids` and
`sequence_ids`; it does not change counting or attribution. Existing JSON schemas
remain stable. Reports should read producer measurements rather than reconstruct
support from assembly voters or count serialized rows.

## What each measurement means

| Native path field | Meaning |
|---|---|
| `rna_support.reads`, `supporting_segments` | Distinct original sequenced segments contributing an exact placed overlap |
| `rna_support.fragments`, `supporting_fragments` | Distinct original templates among those contributors; mates count together |
| `coverage` | Original-segment depth at each query base; overlapping mates can contribute two segments |
| `compatible_fragments` | Templates whose observed segments are compatible with this candidate |
| `unique_fragments` | Compatible templates assigned to only this retained candidate, with no unresolved outside placement |
| `whole_path_segments` | Original segments matching the complete path |
| `coding[].whole_cds_segments` | Original segments matching that complete coding interval |
| `coding[].all_bases_observed` | Every coding base is covered; does not require a single spanning segment |

`rna_support` uses the same read/fragment units as small-variant and SV support.
The *subject* of a count still matters: focal ALT support, a direct junction,
partial path contributors and whole-CDS witnesses are different measurements.
Small-variant `coding_support.fragment_depth` counts fragments, whereas native
path `coverage` counts segments; a report must label that distinction.

Alternative isoforms may share the same compatible fragments. Their counts are
not additive. An ambiguous assignment remains an equivalence class instead of
being split into invented unique counts. Paired ends do not fill an unsequenced
insert, and two tiled reads do not become one end-to-end witness. These are RNA
observations, not counts of independent molecules, protein expression or antigen
presentation.

The shared `rna_evidence.evidence_set` and `union_rna_support` APIs retain scoped
original identities when combining evidence. Different original libraries remain
separate; reused processing views are not independent replication. Unmeasured,
unavailable and measured-zero outcomes remain distinct.

## Coordinates, placement choices and cache ownership

- A `Position` is `(contig, zero_based_base, strand)`. `None` means an observed
  base has no resolved placement; it does not mean a missing base was supplied.
- `PathPositions` is an immutable placement tuple with cached geometry.
  `position_runs` returns half-open query and genomic intervals; genomic bounds
  are low/high on either strand. `path_breaks` returns the two *flanking base
  indices*, not an interval.
- `PlacementIndex` routes shared placed bases. A routed overlap still requires
  `exact_overlap_agrees`; routing alone is not evidence of sequence support.
- `PathSpans` returns a half-open co-observed geometric span. Placed-base alleles
  may differ, so it must not substitute for an exact full-sequence witness.
- One `RnaObservation` describes one original segment and placement choice.
  Reciprocal orientations and alternative placements retain that identity.
  Assembly never lets alternative placements of a segment extend one another as
  independent witnesses.
- `ObservationGeometryCache` shares geometry only. Original bases, quality,
  metadata and fragment membership remain attached to each observation.
- `original_library_context` closes resources it creates. A supplied open
  `RnaLibraryIndex` remains owned by its caller and is validated for source,
  reference and eligibility compatibility before reuse.

## Related contracts

- [Read collection and RNA support records](read-collection.md)
- [Coding-window coverage and protein hypotheses](protein-hypotheses.md)
- [Literal-allele RNA paths](variant-rna-paths.md)
- [SV reconstruction and junction evidence](sv-rna.md)
- [CTA and other antigen isoforms](antigen-isoforms.md)
- [Viral paths and coding references](viral-sequences.md)

Shared helpers have public names in their owning modules. Package-level tests
reject cross-module private imports and dependencies from this shared engine
back into antigen adapters. Existing public workflow entry points remain
available; previously private helper names are not compatibility APIs.
