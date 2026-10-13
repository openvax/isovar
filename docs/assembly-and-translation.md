# Assembly and translation

A `ProteinSequenceCreator` controls how reads are assembled into coding
sequences, placed in a reading frame and grouped into proteins. The defaults are
shown:

```python
from isovar import run_isovar, ProteinSequenceCreator

protein_sequence_creator = ProteinSequenceCreator(
    # Peptide size K used to score context; the default target length is 2*K-1.
    protein_context_peptide_length=25,
    # None derives the target from the peptide size (49 aa for K=25).
    protein_sequence_length=None,
    # "balanced", "support" or "context"; see protein context selection.
    protein_sequence_preference="balanced",
    # Balanced mode keeps candidates with at least this fraction of the best
    # candidate's compatible read support.
    min_protein_sequence_support_fraction=0.85,
    # Minimum number of reads covering each base of the coding sequence.
    min_variant_sequence_coverage=2,
    # Bases of reference transcript the cDNA must match before the variant
    # to establish a reading frame.
    min_transcript_prefix_length=10,
    # Mismatches allowed between the cDNA and the reference transcript.
    max_transcript_mismatches=2,
    # Also count mismatches after the variant toward max_transcript_mismatches.
    count_mismatches_after_variant=False,
    # Ranked protein sequences kept per variant; 0 keeps all.
    max_protein_sequences_per_variant=1,
    # Assemble overlapping reads; if False each sequence comes from one read.
    variant_sequence_assembly=True,
    # Minimum overlap, in nucleotides, before two reads are combined.
    min_assembly_overlap_size=30)

isovar_results = run_isovar(
    variants="cancer-mutations.vcf",
    alignment_file="tumor-rna.bam",
    protein_sequence_creator=protein_sequence_creator)
```

Most of these settings are also options of `isovar run` and the protein
commands, such as `--min-variant-sequence-coverage`; `isovar run --help` lists
them.

## How reads become a translated sequence

- Each read keeps its aligned exon blocks and splice junctions, so it is
  compatible with a particular set of annotated transcripts.
- Overlapping reads are assembled only along transcripts they share, and the
  cDNA is translated only in those transcripts' frames.
- A read that ends before an isoform-distinguishing junction stays ambiguous:
  it supports every compatible isoform without being counted twice.
- The frame is carried through each read's own alignment, so an upstream indel
  in the RNA shifts it ([details](aligned-reading-frame.md)).

## How much protein context

For peptide size K, Isovar targets 2K−1 residues (15-mers → 29 aa, 25-mers →
49 aa), enough for every K-mer overlapping a centered single-residue change.
The default `balanced` preference maximizes mutation-overlapping peptide
windows among candidates with at least 85% of the best candidate's compatible
read support. `support` ranks by read support first; `context` ignores the
support budget.

Actual context depends on RNA coverage, and no reference sequence fills missing
RNA. [Protein context selection](protein-selection.md) gives the exact rules.

## Other variants in the assembled RNA

An assembled cDNA can differ from the reference transcript at positions besides
the variant: a nearby germline SNP, another somatic mutation, or an RNA or
sequencing difference. Each `ProteinSequence` groups these edits by origin:

| Property | Edits |
|---|---|
| `known_somatic_transcript_edits` | This variant, plus other input variants in the same RNA (co-somatic) |
| `known_germline_transcript_edits` | Variants from a matched normal, given as `run_isovar(germline_variants=...)` or `--germline-vcf` |
| `unexplained_transcript_edits` | Everything else |

An edit counts as a supplied variant only when it is that variant on the
transcript, including an indel shifted within a repeat. A run of adjacent
changes is split into the supplied variants it contains. Origin is never guessed
from allele fraction.

For the top protein, the `isovar run` table reports
`num_co_somatic_variants_in_top_protein_sequence`,
`num_germline_variants_in_top_protein_sequence` and
`num_unexplained_edits_in_top_protein_sequence`. These work as filters too, for
example `max_num_unexplained_edits_in_top_protein_sequence`.

When callers disagree about the variant itself,
[`isovar allele-interpretations`](allele-interpretations.md) compares their
versions with the reads first.
