# Nuclear versus mitochondrial RNA origin

`assess_rna_origin` and `isovar rna-origin` compare the context of an observed
allele against explicit mitochondrial and nuclear RNA paths. A nuclear path
can represent a nuclear mitochondrial insertion (NuMT), including its nuclear
flanks or a spliced transcript. The result is a **panel-relative heuristic**:
`nuclear_favored` or `mitochondrial_favored` means the observed sequence fits
that origin better among the supplied competitors. It is not proof of origin.

This opt-in assessment uses the shared original-library index, original segment
observations, supplementary-alignment validation and scoped RNA evidence sets.
It does not change translation tables, candidate admission or counts in other
exports. Standard translation still defaults to NCBI table 1; mitochondrial
code selection remains the responsibility of the reference-context resolver.

## Supply competing paths

Each `RNAOriginHypothesis` requires:

- A unique `identifier` and an explicit `origin`: `nuclear` or `mitochondrial`.
- An ACGT `sequence` in path orientation, with enough context for the reads.
- One `positions` entry per base: `(contig, zero_based_position, strand)`,
  or `None` for inserted bases. Contig names must match the alignment header.
- A nonempty, half-open `target_interval` in that sequence. An indel or junction
  target includes placed anchors; a wholly unplaced insertion is rejected.

Supply at least one hypothesis of each origin. Include plausible alternative
NuMT loci, nuclear transcript isoforms, mitochondrial haplotypes and known
patient-specific insertions where available. All hypotheses must carry the
**same target sequence**, allowing reverse complement. For a variant, place
the nominated allele in each homologous target; the focal allele is excluded
from origin scoring. This prevents the variant itself from acting as an origin
marker. An unmatched target is not an origin call.

Paths are RNA sequences: join exons in transcript order and retain their actual
genomic positions. Do not compare a spliced read to an unspliced nuclear genomic
window and interpret its poor fit as evidence for mitochondrial origin. Include
sufficient flanks on both sides; the comparison refuses to favor an origin
when any supplied competitor is too short for the observed context.

The following constructs a panel from two supplied, positive-strand, contiguous
201-base windows. Both text files contain ACGT only, carry the nominated allele
at offset 100, and use the same reference build as the BAM. Substitute the
actual contigs, coordinates and competing sequences; these example coordinates
do not identify a real NuMT.

```python
import json
from pathlib import Path
import pysam
from isovar import RNAOriginHypothesis, assess_rna_origin

def window(identifier, origin, filename, contig, start):
    sequence = Path(filename).read_text().strip().upper()
    return RNAOriginHypothesis(
        identifier=identifier, origin=origin, sequence=sequence,
        positions=[(contig, start + i, "+") for i in range(len(sequence))],
        target_interval=(100, 101))

panel = [
    window("mt-allele", "mitochondrial", "mt-window.txt", "chrM", 1000),
    window("numt-allele", "nuclear", "numt-window.txt", "chr1", 100000),
]
with pysam.AlignmentFile("tumor-rna.bam") as bam:
    result = assess_rna_origin(bam, panel, sample_id="tumor-1", source="rna-library-1")
Path("origin.json").write_text(json.dumps(result, indent=2))
Path("origin-panel.json").write_text(json.dumps({
    "schema": "isovar.rna_origin_panel.v1",
    "hypotheses": [h.to_dict() for h in panel],
}, indent=2))
```

The equivalent CLI uses that panel and writes a new output file:

```sh
isovar rna-origin --bam tumor-rna.bam --hypotheses origin-panel.json \
    --sample-id tumor-1 --source-id rna-library-1 --output origin-cli.json
```

An optional panel `reference_name` records the declared reference identity.
`--genome-fasta` supplies a CRAM decoding reference. Repeated BAM assessments can
reuse `--library-cache`; Python callers can pass an opened `RnaLibraryIndex`
through `library_index`. See [original RNA libraries](rna-library-index.md).

## How evidence is compared

1. Select actual original segments whose observed target bases and placements
   match a nominated path. Retain secondary placements and MAPQ-zero reads by
   default. MAPQ 255 is unavailable mapping quality, not an origin vote.
2. Require target bases at Q20. Mask low-quality or unresolved flank bases.
   Missing base qualities or fewer than 30 reliable flank bases leave the
   observation unassessed.
3. Compare the entire eligible observed flanks against each competitor,
   anchored at the target. Edit distance permits substitutions and indels;
   only distal reference overhangs are free. The best origin must have no more
   than 5% edits relative to reliable context bases and beat the other origin
   by at least three edits **plus one per masked context base**. The extra margin
   prevents unknown bases from providing apparent certainty. These thresholds
   are engineering defaults, not a calibrated likelihood or error rate.
4. Combine alternatives and target-bearing mates conservatively. Opposing
   favored origins within one fragment are conflicting. A favored placement
   plus an uncertain alternative remains ambiguous. Repeated SAM records and
   alignments never create additional read or fragment votes.

`RNAOriginPolicy` exposes `min_base_quality`, `min_context_bases`,
`max_edit_fraction` and `min_edit_margin`. CLI flags use the same names with
hyphens. Every result records the effective policy, comparisons and reasons.

Actual linked supplementary records can connect the target to distinguishing
nuclear sequence. An SA tag without its actual compatible partner cannot do
so. In this first version, a disjoint mate does not establish allele-to-context
phase: only segments that themselves contain the target contribute. No assembly
across independent fragments is used to invent origin linkage.

## Read the output

`isovar.rna_origin.v1` includes the panel and its content identifier, input and
read-policy provenance, per-fragment/segment/observation assessments, edit
distances, and source-scoped evidence-set identifiers. It does not expose raw
read names. `rna_support` counts eligible target-bearing reads and fragments;
`support_by_status` partitions those counts into disjoint fragment categories:

| Status | Meaning |
|---|---|
| `nuclear_favored` | Context favors a supplied nuclear path |
| `mitochondrial_favored` | Context favors a supplied mitochondrial path |
| `ambiguous` | Tied or insufficiently separated origins, or uncertain alternatives |
| `conflicting` | One fragment has opposing favored origins |
| `unassessed` | Inadequate quality/context, incomplete panel context, or poor fit |

At sample level, `mixed_origins` preserves separately favored nuclear and
mitochondrial fragments. `no_target_support` means no eligible target-bearing
segment was found. Neither state implies a DNA genotype. Counts cover the
declared input, which may be a regional BAM, and are not whole-library abundance,
DNA heteroplasmy or an estimated NuMT contamination fraction. Libraries are
assessed separately; read groups and mate identities retain their shared-engine
meaning. No UMI deduplication is inferred.

## Limits and scientific basis

The caller supplies the panel; this command does not discover genome-wide or
patient-specific NuMTs. A missing near-identical NuMT can still make mitochondrial
origin look favorable. Short shared sequence may be fundamentally unidentifiable.
RNA library chimeras can also mimic molecular linkage. Independent DNA junction
or long-read evidence is valuable when an origin distinction matters.

Nuclear–mitochondrial junctions and long-read validation are supported by
[Wei et al., 2022](https://pmc.ncbi.nlm.nih.gov/articles/PMC9630118/).
Allele-dependent cross-mapping motivates masking the focal allele and comparing
longer context ([Maude et al., 2019](https://doi.org/10.3389/fcell.2019.00201)).
Non-reference NuMTs and copy-number-dependent false positives limit what a
reference panel can exclude ([Laricchia et al., 2022](https://pmc.ncbi.nlm.nih.gov/articles/PMC8896463/)).
The thresholds above have synthetic regression coverage, not biological
calibration. Native-DNA methylation is outside this RNA/cDNA assessment.
