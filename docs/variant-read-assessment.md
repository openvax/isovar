# RNA assessment before protein selection

An ALT allele can have assembled RNA without an assignable coding frame.
`ProteinSequenceCreator.assess_variant_reads` preserves both results. Intronic,
noncoding and coding literals use the same assembly policy; absence of a coding
annotation no longer skips assembly.

```python
from isovar import ProteinSequenceCreator, ReadCollector

creator = ProteinSequenceCreator(max_protein_sequences_per_variant=0)
evidence = ReadCollector().read_evidence_for_variant(variant, bam)
assessment = creator.assess_variant_reads(variant, evidence.alt_reads)
proteins = creator.protein_sequences_from_assessment(assessment)

print(assessment.assembly.state, assessment.frame_state)
for sequence in assessment.assembly.unmatched_sequences:
    print(sequence.sequence, len(sequence.reads))  # contributors, not full-span support
```

`VariantReadAssessment` and `VariantRNAAssembly` are frozen containers with typed
fields. They reference the existing `VariantSequence`, `AlleleRead`,
`ReferenceContext` and `Translation` objects. They introduce no mutable
last-result state on the creator. Ranking an assessment does not repeat assembly
or remove its alternative paths. Its `translations` precede protein selection;
synonymous/reference-identical contexts remain subject to the existing separate
protein-comparison and mutant-count rules.

## Interpret the states

| Field | State | Meaning |
|---|---|---|
| `assembly.state` | `assembled` | At least one local ALT RNA candidate passed the configured assembly and coverage policy. |
| | `no_alt_reads` | No ALT observations were supplied. |
| | `no_retained_rna` | ALT observations were supplied, but none produced a retained candidate under these budgets. |
| `frame_state` | `rna_frame_candidates` | At least one retained RNA context matched a supplied reference frame. |
| | `no_coding_reference_context` | Ordinary-CDS lookup supplied no usable frame at this allele. |
| | `no_reference_path_compatible_reads` | Observed alignments conflict with every supplied coding transcript path. |
| | `no_retained_compatible_rna` | No compatible RNA candidate survived assembly and coverage. |
| | `no_retained_rna_frame_match` | Compatible RNA survived, but no sequence/frame match passed the matching policy. |

Assembly and frame states are independent: a missing frame is not missing RNA.
The `transcript_id_whitelist` restricts frame hypotheses without suppressing
unmatched RNA. RNA compatible with an annotated splice path and RNA incompatible
with all supplied paths are assembled separately, with identical coverage and
context budgets. The latter is never forced into that annotated path. An
explicit `ReferenceContext` with no transcript metadata may still supply a
conditional sequence frame; its empty transcript IDs make the absent path
constraint explicit.

Use `assembly_mode="read_contexts"` to retain bounded original read contexts
without overlap extension. The cohort and evidence-report commands use this
same route at their high-depth cutoff. Both routes now apply the creator's
protein selection and count cap consistently. The default `"configured"`
uses the creator's existing assembly and context-length settings.

## API, CLI and stored results

`run_isovar` attaches the assessment to `IsovarResult.rna_assessment`. Cloning and
pickling preserve it. The [protein-hypotheses JSON](protein-hypotheses.md) stores
an event-level `rna_assessment`, including local nucleotide sequences,
zero-based half-open variant intervals, per-path frame states, compatible and
frame transcript IDs, and contributor evidence-set IDs. These share the event's
sample/library identity and do not create independent read support. Protein TSV
rows remain protein/translation rows; unframed RNA is available in JSON.

Older results and hand-built results without an assessment export `null`, which
means unassessed. It does not mean that assembly failed. A custom creator that
overrides the older list-only translation/ranking hooks retains control of its
protein output and likewise has no inferred RNA assessment. To customize the
complete flow, override `assess_variant_reads` and/or
`protein_sequences_from_assessment` instead. The existing
`translate_variant_reads` and `sorted_protein_sequences_for_variant` list APIs
remain available.

[Supplied ORFdata frames](variant-rna-orfs.md) use the same assembly helper and
offer the same `assembly_mode`. Their model-specific outcomes remain separate
from ordinary-CDS assessment. An absent ncORF catalogue is not a reason to omit
the ordinary assessment's unframed RNA. Reference contigs are resolved once with
Varcode: exact custom/viral/decoy names win, ambiguous/missing names remain
explicit, and mitochondrial aliases do not rewrite arbitrary substrings.

## Scope and biological interpretation

This is local, allele-centred assembly, not a whole-transcript search. Each
candidate retains its original witnesses and placement alternatives. A
contributor count is not a count of reads spanning an entire RNA or protein.
Where a coding candidate exists, use its existing interval-specific coverage,
single-read and paired-template support fields.

Retained-intron RNA can yield peptides in some settings ([Smart et al., 2018](https://doi.org/10.1038/nbt.4239)).
Depth at an intronic allele alone does not establish productive splicing,
initiation, translation, or a mutation-induced isoform. An
[Ensembl retained-intron model](https://www.ensembl.org/info/genome/genebuild/biotypes.html)
does not by itself supply a CDS. Recovering original reads that link an intronic
allele to a coding exon or cryptic junction, comparing their alternative
placements, and assigning a justified partial/complete ORF remain separate
evidence requirements tracked in #583. Reference sequence is never inserted as
observed RNA.
