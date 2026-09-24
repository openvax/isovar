# Supplied fusion RNA

`isovar fusion` checks a fusion transcript that another tool has already
assembled. It asks: *do the RNA reads actually show this sequence and its
junction, and does an annotated reading frame carry across the junction into a
new protein?* It does not discover fusions, fill gaps with reference sequence,
pick the longest ORF, or assemble a transcript. A gene pair, a symbolic VCF
allele or a DNA breakpoint alone is not a fusion transcript.

To find RNA paths around an SV call from a BAM, without a transcript in hand,
use [`isovar sv-rna`](sv-rna.md). Both commands return the same
[RNA path format](sv-rna.md#the-rna-path-format). Vaxrank reads the result
with `fusion_antigens_from_isovar`; see
[library responsibilities](library-responsibilities.md) for how this fits with
Varcode and Vaxrank.

## Quick start

```sh
isovar fusion --input fusion.json --output result.json
# Optional figures in a new UTC-stamped directory:
isovar fusion --input fusion.json --output result.json --plot-dir figures/fusions
```

From Python:

```python
from isovar import reconstruct_fusion
from isovar.fusion import fusion_from_dict

fusion, references, reads = fusion_from_dict(data)
result = reconstruct_fusion(fusion, references, reads)
```

## Input

The input JSON has these keys; unknown keys are rejected:

| Key | Required | Content |
|---|---|---|
| `fusion` | yes | The supplied transcript (`FusionTranscript`) |
| `references` | no | Every candidate transcript model for both partners (`FusionReference`), not just the preferred one |
| `reads` | no | RNA observations of the transcript (`FusionRead`) |
| `reference_names` | no | Transcript ID → display name, used only to label figures |

A minimal `fusion`:

```json
{
  "event_id": "GENE1--GENE2",
  "reference_name": "GRCh38",
  "sequence": "ATGGCTGCTAAACCTGGGTCCCTTTGAATAA",
  "junction_start": 15, "junction_end": 15,
  "donor":    {"contig": "1", "position": 1115, "strand": "+"},
  "acceptor": {"contig": "2", "position": 2003, "strand": "+"},
  "blocks": [
    {"query_start": 0,  "query_end": 15, "contig": "1", "reference_start": 1100, "reference_end": 1115, "strand": "+"},
    {"query_start": 15, "query_end": 31, "contig": "2", "reference_start": 2003, "reference_end": 2019, "strand": "+"}
  ],
  "provenance": {"sample_id": "tumor-1", "method": "assembler", "version": "1.0",
                 "parameters": {}, "source": "rna.bam", "contig_id": "contig-7"}
}
```

- **Coordinates** are 0-based and half-open. Genomic intervals are ascending,
  even on the minus strand. Positions within the transcript (`query_*`, CDS
  offsets) run 5′→3′ along the supplied RNA.
- **`junction_start:junction_end`** is the inserted sequence between the partners;
  it is empty for a direct join.
- **`blocks`** align the two partner sequences to the genome. They must reach the
  junction and agree with the `donor`/`acceptor` breakpoints, which are each
  partner's retained boundary, not a VCF position.
- **`provenance`** must name the sample, method, version, parameters, input source
  and contig ID. It is kept verbatim.

Each **read** gives `sample_id`, `library_id`, `fragment_id`, `read_id`, `source`,
its `sequence`, where it sits on the supplied transcript (`cdna_start`), its own
alignment `blocks`, `source_query_start` in the original read, and `mate_number`
(1 or 2 for paired ends, 0 otherwise). Resolve supplementary pieces into one
observation per mate before submitting, and never union alternative placements.
Copies of one read from different processed files are counted once; conflicting
copies are an error.

**References** give each transcript's exons, spliced sequence and optional CDS.
The CDS end includes the stop codon, and a missing CDS means noncoding or
unknown. Only NCBI genetic code 1 is supported.

## Result

The result (schema `isovar.fusion_rna.v2`) has one path, `paths[0]`, the
supplied transcript, in the same structure as `isovar sv-rna`:

| Field | Meaning |
|---|---|
| `status` | `translated`, `ambiguous`, `unresolved` or `insufficient_support` (below) |
| `reasons` | Why a frame could not be used, per transcript model |
| `paths[0].sequence`, `sequence_sha256`, `blocks` | The supplied RNA and its alignment |
| `paths[0].junctions[0]` | `query_interval` of the inserted bases, flanking genomic positions `left`/`right`, `unplaced_bases`, `relation` (`breakpoint_junction`) and `direct_fragments` |
| `paths[0].frame_status` | As `status`, or `not_assessed` when support was insufficient |
| `paths[0].translations` | Every justified protein hypothesis |
| `paths[0].compatible_transcripts` | Reference models whose exons exactly match each partner |
| `evidence` | `reads`, `fragments`, `direct_fragments` and every observation, marked with whether it spans the junction |
| `parameters` | `peptide_lengths`, `min_fragments` and `genetic_code` |

### Status

- **`insufficient_support`**: fewer than `--min-fragments` (2) distinct fragments
  directly span the junction, or some transcript bases are not covered by any
  read. A supplied sequence alone is never a validated protein.
- **`unresolved`**: the RNA is supported, but no annotated frame can be justified.
  This does not mean no fusion RNA exists.
- **`ambiguous`**: more than one protein hypothesis, or coding and noncoding donor
  models both fit. Downstream tools must not silently pick one.
- **`translated`**: one protein hypothesis across the compatible donor models.
  This is not a uniquely chosen isoform or proof of protein expression. Check
  `complete_5prime` and the frame evidence before use.

### Translations

A frame is transferred only when the donor side matches an annotated coding
transcript exactly and collinearly. Upstream indels, unknown starts, novel
donor splicing and sequence mismatches leave the frame unresolved. The supplied
RNA is then translated through its first stop or its end. Each translation has:

- `translation_start`/`translation_end`, `amino_acids` and `ends_with_stop_codon`;
- `complete_5prime`: whether the annotated start codon is in the RNA. When false,
  `cds_start` is null and the protein assumes the annotated upstream frame;
- `transcript_ids` and `frame_evidence` for every donor model giving this protein;
- `acceptor_frames`: whether each acceptor model is in frame;
- `candidate_peptides`: every window of the requested lengths that crosses a
  junction boundary, with `junction_boundaries_in_cds`. An insertion has two boundaries;
- `downstream_frameshift_peptides`: windows wholly after the junction in a frame
  different from an acceptor model;
- `junction_in_translated_cds` and `trailing_partial_codon_bases`.

Neither peptide list establishes absence from the reference proteome,
presentation or immunogenicity.

## Figures

`--plot-dir DIR` writes white-background PNG and SVG panels and a vector PDF to a
new UTC-stamped directory under `DIR`, never overwriting an earlier run. There
are separate panels for the junction reads, local donor and acceptor
annotation, and each justified protein. `--dpi` sets PNG resolution (default
600, minimum 72). The JSON result is written before plotting, so a plotting
error does not lose it.

## Real RNA examples

The pinned osteosarc examples keep original ONT records and the exact repeated
junction windows, not reference-built transcripts. TPST1–CRCP has a donor 5′-UTR
junction. The FOXO3 and PARD3B windows are intronic relative to the supplied
Ensembl 87 donor models. All of them stay `unresolved` rather than acquiring an
invented coding fusion. The PARD3B window keeps its observed 12-nt insertion.
Counts describe the selected sequence group, not total event abundance. Fusion
results must not be turned into a single-locus `Variant` for Vaxrank.

Sources: [SAM format](https://samtools.github.io/hts-specs/SAMv1.pdf),
[NCBI translation tables](https://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi),
[osteosarc fusion evidence](https://osteosarc.com/fusions/),
[CTAT-LR-fusion](https://github.com/TrinityCTAT/CTAT-LR-fusion/wiki).
