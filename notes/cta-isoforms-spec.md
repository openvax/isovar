# CTA isoforms and a combined alignment reference

Implement the first usable, variant-free part of #465: reconstruct targeted
RNA paths, retain competing annotated isoforms, report original fragment
compatibility, and group complete encoded proteins. A caller supplies the CTA
list; Isovar does not silently define cancer-testis specificity.

Use the existing coordinate-aware RNA observation and exact-overlap assembly
code. Never insert reference bases into an observed path. Distinguish reference
models, partially observed sequence, overlap-assembled hypotheses, and paths
covered by individual reads. Keep unsupported coding models and missing
coverage explicit. Discovery limits must not limit the subsequent support scan.

Abundance is a separate, optional Salmon step against the complete supplied
transcriptome augmented with the reconstructed candidates. Retain uniquely and
ambiguously compatible fragments independently of the estimates, report
within-gene transcript fractions, and sum transcript estimates by identical
complete protein sequence. Do not turn a numerical allocation between
indistinguishable candidates into a resolved dominant isoform.

Recommend GDC GRCh38.d1.vd1 for new combined human/viral alignments. It is an
alignment reference, not an Ensembl release. Use PyEnsembl's existing custom
Genome API for GTF annotation and genomic DNA; retain the assembly, alignment
profile, annotation release and file checksums separately. Validate contig names,
lengths and available sequence checksums. Existing variant CLI defaults and
existing BAM/reference identities remain unchanged. Viral contigs without coding
annotation do not become viral transcript models or evidence of absence.

Deliver Python/CLI entry points, documentation, offline regression cases for
competing isoforms/paralogs, complete/partial CDS, sequence alleles, both strands,
secondary placements and mate/library identity, and an actual Salmon mixture
experiment. Then run lint, the full suite and CI, merge the PR and publish the
versioned release. Broader noisy-read inference, viral coding exceptions and
sample-specific viral reference construction remain tracked in #465.

Primary sources checked before implementation: GDC reference files and viral
sequence list; GENCODE GTF format; Salmon (Patro et al., Nature Methods 2017,
doi:10.1038/nmeth.4197). Human GTF/CDS offsets are converted from one-based
inclusive coordinates to the existing FusionReference interbase convention.
