# RNA reconstruction in supplied coding frames

`isovar.reconstruct_variant_orfs` assesses validated ORFdata `TranscriptORF`
models against mutation-supporting RNA. The models may describe ordinary CDSs,
upstream/downstream/internal ORFs, or ORFs on noncoding transcripts. The API
preserves the real parent transcript and does not require an annotated CDS.

```python
import orfdata
import pysam
from isovar import ReadCollector, reconstruct_variant_orfs

# variant is a Varcode Variant with its pinned PyEnsembl genome.
catalog = orfdata.load_catalog()
record = catalog.get("c17norep133")
transcript = variant.genome.transcript_by_id(record.transcript_id)
references = [record.bind(transcript)]  # validates the entire ORF path/protein
if transcript.is_protein_coding and transcript.complete:
    references.append(orfdata.TranscriptORF.from_annotation(transcript))

with pysam.AlignmentFile("original-rna.bam") as bam:
    reads = ReadCollector().read_evidence_for_variant(variant, bam)
result = reconstruct_variant_orfs(
    variant, reads, references, sample_id="COLO829", source="SRR8615617")
```

The JSON-ready result contains `rna_candidates`, all supplied `models`,
`protein_candidates`, scoped `evidence_sets`, budgets and explicit unresolved
reasons. The caller supplies and pins models; this API does not choose a
catalogue or download reads. Structural events use `reconstruct_sv_rna` and
its assessment/reconciliation APIs.

Each protein candidate contains the observed nucleotide sequence and coding
interval, mutant protein, reference window, focal mutation interval, and a
counterfactual protein with only the focal allele replaced on the same observed
RNA background. `focal_amino_acid_change` therefore separates synonymous focal
alleles from changes caused elsewhere in the RNA. Unclassified nonfocal
differences must not be labelled germline or co-somatic without independent
genotype and phasing evidence.

`rna_support` counts reads/fragments participating in assembly. `coding_support`
adds exact per-base fragment depth, tiling completeness, and the separate set of
fragments with a single read spanning the complete displayed coding interval
(including the stop codon when present). Tiled support is not end-to-end
molecular linkage. Alternative placements retain separate candidates without
counting the same molecule again. The configured protein-count cap does not
discard models or retained local RNA alternatives in this API.

`complete_coding_candidate` means that the RNA includes the supplied reference
start and an observed in-frame stop. `initiation_observed` and
`translation_observed` remain false. A changed reference start is unresolved;
the API does not guess a new initiation site. Unsupported flanks are never
filled from reference. Partial references and unmatched RNA/frames remain
explicit. This is local allele-centred assembly, not unrestricted cryptic-exon
or whole-transcript discovery.

The [cohort replay](../reports/hcc1395_colo829/ncorf_rna_audit/README.md)
demonstrates default coverage and separately labelled single-read sensitivity
passes on the same original molecules. Never add support across those passes.
Reference ORF evidence follows [Mudge et al. 2022](https://doi.org/10.1038/s41587-022-01369-0);
it does not establish mutant translation or antigen presentation in a sample.
