# Protein hypothesis export

`isovar protein-hypotheses` exports every mutant protein Isovar reconstructs
for each variant, the nucleotide translations behind each protein, and the RNA
reads that support them. It is meant for tools that combine or re-rank
candidates from several sources, such as Topiary and Vaxrank. It does not pick
the top protein, and it adds nothing Isovar did not observe.

The other table commands report one row per protein or translation. This format
also keeps:

- **every alternative**, in Isovar's ranked order, including synonymous cDNAs
  that give the same protein;
- **which hypotheses are shorter windows of another**, without merging them;
- **which reads support what**, as hashed read identities, so that evidence from
  different hypotheses, runs or exports can be combined without counting a read
  twice.

## Quick start

```sh
isovar protein-hypotheses --vcf variants.vcf --bam tumor-rna.bam \
    --sample-id tumor-1 --output tumor-1.hypotheses.json
```

This writes `tumor-1.hypotheses.json` and `tumor-1.hypotheses.tsv`, with one TSV
row per translation. The command takes the same read, assembly, filter and
`--germline-vcf` options as `isovar run`. It keeps every protein by default
(`--max-protein-sequences-per-variant 0`). `--source` names the read set for the
evidence IDs and defaults to the `--bam` path.

From Python:

```python
from isovar import (
    ProteinSequenceCreator, export_protein_hypotheses, run_isovar, write_protein_hypotheses)

results = run_isovar(
    "variants.vcf", "tumor-rna.bam",
    protein_sequence_creator=ProteinSequenceCreator(max_protein_sequences_per_variant=0))
export = export_protein_hypotheses(results, sample_id="tumor-1", source="tumor-rna.bam")
write_protein_hypotheses(export, "tumor-1.hypotheses.json")
```

Pass `alignment_header=pysam.AlignmentFile(path).header` to record each read
group's sample and library.

## Add hypotheses to an ORFdata Proteome

The Python adapter accepts native results or the saved v2 JSON. Supply the
species, assembly and annotation explicitly: the Isovar export records the
build, but does not prove which annotation release was used.

```python
import json

from orfdata import ReferenceContext
from isovar import from_protein_hypotheses, to_proteome

reference = ReferenceContext(9606, "human", "GRCh38", "Ensembl", "87")
with open("tumor-1.hypotheses.json") as handle:
    rna = from_protein_hypotheses(json.load(handle), reference=reference)

# The equivalent native path, using the same mapping and earlier results:
native_rna = to_proteome(results, sample_id="tumor-1", source="tumor-rna.bam",
                         reference=reference)
assert native_rna.identity == rna.identity

# Given an ORFdata reference_proteome built for this same context:
joint = reference_proteome.merge(rna)
joint.write_fasta("joint.fasta")
joint.save("joint.proteome.json.gz")
```

ORFdata 0.39.1 or newer is required. The adapter keeps one FASTA entry per
retained Isovar protein by default, including shorter windows and equal amino
acid sequences from different sources. `write_fasta(..., deduplicate=True)`
is an explicit sequence-grouping choice. The versioned snapshot retains the
original event, hypothesis, translation and evidence-set IDs, their links,
partial boundaries, filter results and completeness/cap states. An event with
no retained protein still has an event record. RNA sequence support remains
RNA evidence; it does not establish that a protein was translated. The first
FASTA token identifies origin, tool, gene, amino-acid length, sample and read
set. ORFdata [#126](https://github.com/openvax/orfdata/issues/126) tracks
display of additional source metadata and partial boundaries in the header;
those values are already preserved in the snapshot.

`from_protein_hypotheses` covers the small-variant v2 export. The shared
[`from_rna_result` API and `isovar proteome` CLI](rna-proteome.md) also accept
SV ORF, supplied-fusion and variant-associated ncORF results. A contradictory
reference build or RNA evidence sample/library is rejected. Linked support
counts must agree with the referenced evidence set; changing only the export's
sample/source labels cannot reassign existing RNA evidence.

## Reading the export

The file (schema `isovar.protein_hypotheses.v2`) has one **event** per input
variant, in input order. Intervals are 0-based and half-open.

| Event field | Meaning |
|---|---|
| `event_id`, `variant` | A stable ID and the variant as given (`start` is Varcode's 1-based position) |
| `reference_prediction` | Varcode's reference-plus-edit effect, for comparison only |
| `allele_support` | The support of the `ref`, `alt` and `other` alleles and their `total`, and `cells_with_ref_and_alt` |
| `filters` | Filter outcomes as recorded; they are reported, not applied |
| `phased_variants` | Event IDs of variants phased with this one |
| `edit_attribution` | How many supplied somatic and germline variants the edits were checked against; null if none were |
| `protein_sequence_limit`, `protein_hypotheses_complete` | The protein cap in use, and whether the list below is everything Isovar found. `false` means the cap was reached; null means unknown |
| `protein_hypotheses` | Every protein, in Isovar's ranked order |

Each **protein hypothesis** has:

| Field | Meaning |
|---|---|
| `hypothesis_id`, `isovar_rank` | A stable ID and the position in Isovar's ranking. The ranking is a policy, not a selection |
| `amino_acids`, `protein_sequence_id` | The protein, and an ID shared by every hypothesis with the same sequence |
| `mutation_interval`, `mutant_amino_acids` | Codons marked by the focal nucleotide edit; they can encode reference-identical residues. A codon-aligned deletion has an empty interval |
| `n_terminus` | `annotated_start_codon` when translation begins at the annotated start, `local_window` when the protein is only the translated window around the variant, `mixed` when its translations differ |
| `c_terminus` | `stop_codon`, or `no_stop_in_context` when coverage or the requested length ended first |
| `representative`, `covered_by` | Whether no other hypothesis contains this one. `covered_by` lists the longer hypotheses it is an exact window of, with its interval in each |
| `transcript_ids`, `gene_names` | Annotated transcripts that give the reading frame |
| `rna_support` | Reads assembled into this protein |
| `translations` | Every cDNA and reading frame giving this protein |

A shorter hypothesis stays in the list even when it is `covered_by` a longer
one. It can have different reads, and its counts are never added to the longer
one's. Filter on `representative` to see one entry per distinct sequence.

Each **translation** has its `nucleotide_sequence` and ID, the
`translated_interval` of the cDNA, the `variant_cdna_interval`, whether it
`starts_at_annotated_start_codon`, the reference context (strand, transcripts,
start codon and 5′ UTR flags), mismatches against the reference, and
`observed_edits`. Observed edits are the differences from each transcript, each
with an `origin`:

- `nominated_variant`: the event's own variant;
- `co_somatic_variant`: another input variant in the same RNA;
- `germline_variant`: a variant from the matched normal (`--germline-vcf`);
- `unexplained`: anything else.

Edits from a supplied variant give its `source_event_id`. Two synonymous cDNAs
are two translations of one protein, each with its own reads.

## RNA support and evidence sets

Every support is the [RNA support record](read-collection.md#the-rna-support-record):
`reads`, `fragments`, and, with `--cell-umi-labels`, `umis`, `cells` and their
completeness. The alternate allele, each protein and each translation also have
an `evidence_set_id`; the reference, other and total allele supports give counts
only, since their reads are rarely combined and dominate at deep loci. The
export's `evidence_sets` table stores each distinct set of reads once:

```python
support = export["events"][0]["protein_hypotheses"][0]["rna_support"]
reads = export["evidence_sets"][support["evidence_set_id"]]
# reads["read_ids"], reads["fragment_ids"], reads["evidence_scope"]
```

A read's ID is the SHA-256 of its read group, name and mate flags, together with
the **evidence scope** `[sample_id, source]`. No read names are exported. Within
one scope, the same read gets the same ID in every export, including
`isovar sv-rna` ORF exports. So to count the evidence behind several hypotheses
or exports, combine their read sets rather than adding their counts:

```python
from isovar import union_rna_support

combined = union_rna_support(
    [export["evidence_sets"][p["rna_support"]["evidence_set_id"]]
     for p in export["events"][0]["protein_hypotheses"]])
combined["reads"], combined["fragments"]
```

`union_rna_support` refuses read sets from different scopes, because their
overlap cannot be known. Reads reprocessed under another `source`, for example,
may be the same molecules. An `evidence_set_id` is null only when some reads
were built without alignment identities; their counts are then given but cannot
be combined. No reads is the empty set, which has an ID.

With `--cell-umi-labels` (Python: `cell_umi_alignment_file=`), the `umis` and
`cells` of every support are counted, and `allele_support` gives
`cells_with_ref_and_alt`; see [cell/UMI evidence](cell-umi-evidence.md#small-variants). Read and fragment
counts are not molecule counts, abundance or TPM. The
`read_groups` table gives each read group's sample, library and platform from
the BAM header, or null when unknown.

## Coding coverage and sequence comparisons

The existing `rna_support` and `ProteinSequence.num_supporting_fragments` count
**contributors**. A contributor need not span the coding window. Each
translation's `coding_support` measures the exact `coding_interval` on that
nucleotide path, including an observed stop codon:

| Field | Measurement |
|---|---|
| `fragment_depth`, `minimum_fragment_depth`, `fully_tiled` | Distinct exact original fragments at each coding base; tiling does not require one spanning observation |
| `allele_observation_span_support` | Exact full-interval observations, possibly made from overlapping merged mates |
| `single_read_span_support` | Spanning observations made from one original read; conservatively excludes merged observations, without recounting their constituent mates |
| `paired_template_span_support` | Spanning observations made from two overlapping mates; unsequenced inserts never supply bases |

Each support record has scoped witness IDs. Subsets can share a fragment and
are **not additive**. Counts on different nucleotide paths must not be added
to claim support for one assembly. A complete coding candidate additionally
needs the annotated start and observed stop; local window coverage alone does
not establish complete-sequence support or translation.

`reference_comparison` retains two independent comparisons:

- `annotation_comparison`: observed protein versus the matched annotation frame;
- `focal_allele_comparison`: observed protein versus replacing only the focal
  allele on the same observed RNA background.

Both return `status` (`reference_identical`, `amino_acid_divergent`, or
`unresolved`), zero-based substitution/insertion/deletion positions,
`changed_residue_count`, unknown alignment columns, and `stop_changed`.
The residue count includes substitutions, inserted residues and deleted
reference residues. A premature or shifted stop is tracked separately.
Unavailable reference residues leave the count null, even when another
verified difference establishes divergence. An unavailable comparator codon
cannot establish a stop gain. Comparisons are anchored to the nucleotide
edit, not a similarity alignment of repeated amino acids.

A synonymous focal allele can coexist with a changed RNA background. It is
then reference-divergent but not a focal changed protein; unrelated residues
remain `unexplained` unless supplied genotypes and exact edit replay establish
their origin. The marked `mutation_interval` remains available independently
of these comparisons.

Report classification requires compatible models to agree. Conflicting or
missing comparisons remain unresolved. One path is selected for display by
focal change, exact span support, minimum depth and stable identity. Its
`translation_id` and interval identify precisely which coverage is displayed;
the full export retains every alternative. TSV rows expose these comparisons
and contributor, span, paired-template and single-read counts separately.

## What the export does not do

- It never chooses one protein, applies filters or drops a variant without RNA
  support. Those remain later policies.
- It only contains what Isovar kept. With a cap above 0, check
  `protein_hypotheses_complete`.
- A `local_window` protein is not a complete ORF. Neither `n_terminus` nor
  `c_terminus` establishes translation, presentation or tumor specificity.

## Partial-read support scores

`--partial-read-support` adds catalog-relative fractional fragment scores and a
separate score rank to JSON and TSV. It evaluates original partial reads against
all exported candidates while preserving raw counts and Isovar ordering. With
`--cell-reconstruction`, the same score is also summarized per cell.
See [partial-read support](partial-read-support.md) for counting, ambiguity,
quality categories and comparison settings.
