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

## Reading the export

The file (schema `isovar.protein_hypotheses.v1`) has one **event** per input
variant, in input order. Intervals are 0-based and half-open.

| Event field | Meaning |
|---|---|
| `event_id`, `variant` | A stable ID and the variant as given (`start` is Varcode's 1-based position) |
| `reference_prediction` | Varcode's reference-plus-edit effect, for comparison only |
| `allele_support` | Alt reads with their evidence set; ref, other and total as counts |
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
| `mutation_interval`, `mutant_amino_acids` | Where the change is in the protein; a deletion has an empty interval |
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

Every `rna_support` gives `segments` (sequenced segments: a mate or a single
read) and `fragments` (templates: both mates count once), with an
`evidence_set_id`. The export's `evidence_sets` table stores each distinct set
of reads once:

```python
support = export["events"][0]["protein_hypotheses"][0]["rna_support"]
reads = export["evidence_sets"][support["evidence_set_id"]]
# reads["segment_ids"], reads["fragment_ids"], reads["evidence_scope"]
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
combined["segments"], combined["fragments"]
```

`union_rna_support` refuses read sets from different scopes, because their
overlap cannot be known. Reads reprocessed under another `source`, for example,
may be the same molecules. An `evidence_set_id` is null when some reads were
built without alignment identities; their counts are then given but cannot be
combined.

Read and fragment counts are not molecule counts, abundance or TPM. The
`read_groups` table gives each read group's sample, library and platform from
the BAM header, or null when unknown.

## What the export does not do

- It never chooses one protein, applies filters or drops a variant without RNA
  support. Those remain later policies.
- It only contains what Isovar kept. With a cap above 0, check
  `protein_hypotheses_complete`.
- A `local_window` protein is not a complete ORF. Neither `n_terminus` nor
  `c_terminus` establishes translation, presentation or tumor specificity.
