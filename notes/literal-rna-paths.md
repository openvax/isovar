# Original RNA paths around a literal allele (#583)

Add a library/API and CLI consumer for intronic, coding and noncoding literal
variants. Use the existing `rna_reconstruction` observation/extension/support
engine, `ReadCollector` allele normalization, source-scoped original-template
index and typed transcript-path frame assessment. Do not nominate a fake SV or
implement another report-only reconstruction engine.

## Inputs and one flow

An input names the literal variant, assembly, sample/library, explicitly versioned
transcript/ORF models, source provenance and search intervals. BAM and all models
must share an explicit coordinate system. Varcode resolves contig vocabulary;
it never translates positions or infers offsets from a local-reference name.
Library reuse/pinning, read-end profiles and missing/ambiguous contig states use
the existing implementations. Keep each culture/library separate.

Recover complete original templates overlapping the allele, including mates and
alternative/supplementary placements. Call alleles anew on those placements with
the shared collector. Original ALT/REF labels from a different alignment cannot
be inherited as current phase. Convert the collector's unchanged-query allele
footprints to the existing oriented RNA observation coordinates. Allele-bearing
observations seed discovery; the existing bounded regional observations extend
them through exact coordinate-consistent overlaps. ALT and REF controls use the
same budgets. Preserve all placement alternatives and report reached budgets.

Extend the common discovery interface with explicit seed observations rather
than duplicating assembly. The existing variant-free behavior remains the
default. Carry parent segment/placement identity through both orientations and
technical barriers. No reference-filled bases or synthetic read support.

## Output evidence

For every retained path, keep the sequence/placements, original allele seeds,
focal query intervals and permitted observed splice/indel departures. Assess
every supplied model, including ordinary and noncoding competitors. Record
conditional CDS-anchor and retention readings, allele coverage versus upstream
stop, incomplete starts/ends, and unresolved alternatives. A covered interval
alone is not a changed protein: preserve reference/provenance comparison rules.

Recount original input observations using the common support engine. Separate
contributors, per-base segment coverage, single-read entire-path support,
paired-template coverage, and direct allele-to-junction/CDS-anchor linkage.
Overlap-chain continuity is an assembled hypothesis, not end-to-end molecular
phase. Mates may link features without filling an unsequenced gap. Secondary
placements from one segment must never combine to manufacture phase or extra
independent support. Expose missing/unreported placement evidence.

CLI and cohort integration use the same input decoder and library function.
Store RNA/frame outcomes even when no protein can be assigned. Reports consume
these results, show precise source/recovery scope and evidence budgets, and do
not present an unassessed mode as biological absence.

## DDX17 and verification

Replay the already recovered 27,309 SRR8615617 original pairs and their pinned
splice-aware placements. Their STAR BAM is relative to a 50-kb GRCh38 slice:
any derived whole-genome-coordinate input must verify the local FASTA against
the pinned target slice, translate every applicable coordinate-bearing field,
retain original records and a hash/projection receipt, and preserve local-only
alignment uniqueness. No silent alias-based liftover.

Use fresh ALT/REF calls; compare the prior 221-fragment selected-site result
without assuming the newly recovered data must reproduce those labels. The
existing investigation's ordinary/NMD/cryptic junctions and intronic coverage
are regression controls, not preassigned conclusions. Stop-before-deletion
readings remain reference/path predictions unless an observed allele-linked
path justifies a mutant protein. No negative whole-library or translation claim.

Self-review: plus/minus SNV/indel/repeat-normalization cases, retained/cryptic
paths, true exon-linked ALT and REF controls, upstream stops, incompatible
placements/mates, technical joins, contig ambiguity, budget exhaustion and
CLI/API equivalence. Preserve existing variant-free/SV snapshots. Run lint,
focused tests, full suite and exact-head CI; merge and deploy in series order.

Primary sources: MiSplice https://doi.org/10.1038/s41467-020-19307-6;
retained-intron peptides https://doi.org/10.1038/nbt.4239;
SAM coordinates/tags https://samtools.github.io/hts-specs/.
