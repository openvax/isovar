# Internal design

![Isovar's data flow, from variants and RNA reads to protein sequences](https://raw.githubusercontent.com/openvax/isovar/master/isovar_design.png)

Isovar's inputs are one or more somatic variant call (VCF) files and a BAM file
of aligned tumor RNA reads. These objects carry the information from one stage
to the next:

1. [`LocusRead`](../isovar/locus_read.py): a read overlapping a variant locus.
   It keeps the quality and alignment details used for filtering (for example
   MAPQ > 0), which later stages drop.
2. [`AlleleRead`](../isovar/allele_read.py): a filtered read reduced to its
   cDNA sequence before, at and after the variant locus.
3. [`ReadEvidence`](../isovar/read_evidence.py): the `AlleleRead`s at one locus,
   grouped into ref reads (agreeing with the reference), alt reads (agreeing
   with the mutation) and other reads (any other allele). The alt reads are
   used to determine a mutant coding sequence; the ref and other groups are kept
   for filtering.
4. [`VariantSequence`](../isovar/variant_sequence.py): a longer candidate coding
   sequence that `VariantSequenceCreator` assembles from overlapping alt reads,
   together with the reads used to build it.
5. [`ReferenceContext`](../isovar/reference_context.py): the cDNA upstream of
   the variant and which of the {0, +1, +2} reading frames it is translated in.
   Isovar collapses the Ensembl transcripts overlapping the locus into one or
   more of these.
6. [`VariantORF`](../isovar/variant_orf.py) and
   [`Translation`](../isovar/translation.py): a `VariantSequence` placed in the
   reading frame of a `ReferenceContext`, and its translation into a protein
   fragment.
7. [`ProteinSequence`](../isovar/protein_sequence.py): distinct variant
   sequences and reference contexts can give the same translation, so
   `ProteinSequenceCreator` groups equivalent `Translation` objects into one
   `ProteinSequence`. `TranscriptAssemblyEdit` records the transcript-relative
   edits observed in its assemblies.
8. [`IsovarResult`](../isovar/isovar_result.py): a variant and the
   `ProteinSequence` objects associated with it, since one locus can have reads
   that assemble into several incompatible coding sequences. Protein sequences
   are ranked by the configured context/support preference, and the top one is
   `top_protein_sequence`. Allele-support properties such as
   `num_alt_fragments` and `fraction_ref_reads` remain separate from the
   selected protein's compatible support.
