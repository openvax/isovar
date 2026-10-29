# Original RNA paths around literal variants

Use `reconstruct_variant_rna_paths` or `isovar variant-rna-paths` when local ALT
RNA needs a larger exon, splice or retained-intron context. This opt-in workflow
assembles RNA **before** consulting coding annotation. It shares allele calling,
contig resolution, original-library indexing, technical-sequence policy, exact
path extension and frame assessment with the existing library. Missing CDS
models do not discard an observed RNA path.

## One input for API, CLI and cohort

```json
{
  "variant": {"contig": "chr22", "start": 38489856, "ref": "GA", "alt": "G"},
  "reference_name": "GRCh38",
  "sample_id": "COLO829",
  "references": [],
  "regions": [["chr22", 38485000, 38507000]],
  "provenance": {"scope": "Original RNA from the declared library"}
}
```

The variant uses Varcode's one-based literal-allele convention; regions are
zero-based, half-open. `references` contains versioned `FusionReference` objects
serialized as dictionaries: transcript ID, annotation, assembly, contig, strand,
ascending genomic exons, spliced sequence, and optional validated CDS bounds.
Supply coding, noncoding and credible alternative models together. An empty
list permits RNA discovery but leaves protein interpretation unassessed.

```bash
isovar variant-rna-paths --input allele.json --bam original.bam \
  --source SRR8615617 --output paths.json.gz --library-cache rna-indexes
```

`--max-discovery-reads`, `--max-seed-segments`, `--max-paths`, `--min-overlap`,
`--min-anchor-bases` and `--flank` expose the same API parameters. Read filters,
soft clips and read-end profiles use the common RNA arguments. The default
path policy retains secondary placements and MAPQ-zero observations as explicit
competitors, excludes duplicate-marked reads, and excludes soft-clipped bases.
It does not interpret MAPQ-zero RNA as uniquely assigned evidence. BAM files
need no coordinate index; the original-template index supplies the queries.
A persistent cache requires a seekable BAM. SAM/CRAM use a private context;
CRAM may need `--reference-fasta`.

```python
import pysam
from isovar import reconstruct_variant_rna_paths, variant_rna_path_input_from_dict

arguments = variant_rna_path_input_from_dict(input_dictionary)
with pysam.AlignmentFile("original.bam") as bam:
    result = reconstruct_variant_rna_paths(bam, **arguments, source="SRR8615617")
```

For repeated events, pass one open `RnaLibraryIndex` as `library_index`; see
[original-library reuse](rna-library-index.md). The cohort manifest family
`variant_rna_paths` uses the same JSON input and library function:

```json
"variant_rna_paths": [
  {
    "id": "DDX17", "source_id": "colo-library", "sample_id": "COLO829",
    "input": "ddx17-input",
    "parameters": {"max_paths": 128, "min_overlap": 30},
    "read_collection": {"min_mapping_quality": 0}
  }
]
```

Declare the input JSON and source BAM, including their SHA-256 values, in the
manifest's `files`. Each source needs its own library/culture identity and
scope. Cohort events share one prepared index for an identical source,
reference declaration and read policy. Contexts close on success and failure.
Resume and render use the existing scientific-identity and result-digest checks.

## How to interpret the result

`isovar.variant_rna_paths.v1` retains all selected paths and conditional frames.
ALT and REF discovery have equal, independent budgets; their paths can overlap
and are not independent molecules. All eligible original allele observations
and the support recount remain available even when discovery is truncated.
Original read-group, query-name, mate and placement identities are preserved.

`reference_paths` contains separately labelled ordinary model competitors, with
support recounted from the same input. These supplied sequences are not observed
assemblies or mutant proteins. An identical reconstructed path shares its support
ID and gains `reference_transcript_ids`; it is not counted a second time.

| Evidence | Meaning |
|---|---|
| `allele_counts` | Original segments and fragments called ALT, REF or OTHER. Alternative placements can put one parent in more than one state. |
| `focal_alleles` | Actual CIGAR query intervals, repeat-normalized intervals, original records and approved observation keys. |
| `coverage` / `supporting_segments` | Per-base or contributing original segments; neither is whole-protein coverage. |
| `whole_path_segments` | One segment's single compatible geometry covers the entire RNA path. |
| `linkage.single_read_whole_path_fragments` | Distinct templates containing a complete-path single read. |
| `linkage.paired_template_whole_path_fragments` | Two different ends of a template cover the whole path with no unsequenced gap; can overlap the single-read count. |
| `linkage.paired_template_only_whole_path_fragments` | The paired count after excluding templates with a whole-path single read. |
| `linkage.witnesses[].links` | Direct-read versus paired-template links between focal alleles, path junctions and exact frame anchors. A paired link can span an unsequenced gap and is not continuous coverage. |
| `coding[].focal_coverage` | Whether a conditional frame includes the allele, stops before it, or starts downstream. |
| `protein_interpretation` | Missing models, no justified CDS frame, or conditional frames available; never a verified mutant-protein count. |

The path's actual allele must have passed the common allele caller, including
its quality policy. A similar sequence or an SA tag alone cannot create an
allele witness. Alternate placements of one segment cannot be combined into a
single read or a read pair. Different pairwise feature links need not share one
global placement choice and cannot be assembled into end-to-end molecular phase.

Frames transfer from exact coding anchors over observed sequence only. They
can follow an observed exon-to-intron continuation or a path discontinuity;
`assembled_path_discontinuity` is explicitly conditional, with independent
original-read linkage evidence. Noncoding models remain noncoding competitors.
The full reference protein is exported separately. No missing sequence is
filled from annotation, and no longest ORF substitutes for an absent CDS.
An upstream frame that stops before the allele does not encode that allele.

This result deliberately does not classify protein novelty or variant causation
from RNA allele coverage. Reference-identical substitutions, ordinary isoforms,
retention without mutation-induced splicing, and unobserved translation remain
possible. Use the existing [protein comparison rules](protein-hypotheses.md)
when interpreting candidates; conditional frames must not inflate mutant-protein
funnels. The RNA/splicing distinction follows the [MiSplice study](https://www.nature.com/articles/s41467-020-19307-6).

## Source scope and limitations

Contig aliases use Varcode against the actual BAM vocabulary. Missing or
ambiguous names raise `UnresolvedContigError`, not zero RNA support. Aliases
cannot establish assembly equivalence or translate coordinates. In particular,
a BAM aligned to a local reference slice needs an explicitly validated projection
or models and alleles expressed in that local coordinate system. The input
provenance must retain the original source, reference and alignment scope.

Exact reconstruction can fragment noisy long reads or omit a weak branch under
the configured budgets. Read the separate ALT/REF discovery receipts. Missing
alternative alignments, mates and MAPQ uncertainty remain in support witnesses.
A regional BAM cannot establish whole-library absence. RNA sequence and coverage
alone cannot distinguish immature RNA from a translated mature isoform.

## Replaying the retained COLO829 DDX17 evidence

The repository's original DDX17 audit contains recovered original pairs aligned
with STAR to a 50-kb genomic slice. The preparation example validates those
source checksums, verifies the local reference against the supplied GRCh38 FASTA,
projects the explicitly declared coordinates, preserves the SAM payload, and
writes the standard input and cohort manifest. It does no assembly or translation.

```bash
python examples/prepare_ddx17_rna_paths.py \
  --original-audit reports/hcc1395_colo829/ddx17_audit \
  --genome-fasta /path/to/GRCh38_full_analysis_set_plus_decoy_hla.fa \
  --output-dir ddx17-inputs
isovar cohort --manifest ddx17-inputs/manifest.json --output-dir ddx17-report
```

Use a new output directory. Ensembl 87 supplies annotated CDS boundaries; its
local cache is used or the standard annotation download is required. The
manifest preserves the original recovery interval and local-alignment scope.
It cannot establish whole-genome uniqueness or whole-library splice absence.
Read `allele_counts` separately from each state's `seed_segments`: ambiguous
flanking bases can prevent assembly while leaving the locus allele callable.
