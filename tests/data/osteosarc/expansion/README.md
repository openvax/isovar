# Expanded osteosarc RNA/protein audit (#218)

This is a read-evidence and local protein-reconstruction benchmark, **not** a
somatic-calling, immunogenicity, vaccine-construct or clinical efficacy truth
set. Vaccine inclusion chooses the cohort; it does not establish correctness.
The public sequencing data are CC0 according to the
[AWS registry](https://registry.opendata.aws/sid-osteosarc/).

## Results at a glance

The [#227 performance follow-up](performance/README.md) separately records
completed reruns of all ten previously interrupted mitochondrial modes and
unchanged results for 215 nuclear control rows. Its audit budget is explicit;
the original matrix below remains the historical 1.8.1 snapshot.

The [variant summary](audit/SUMMARY.md) and [source summary](audit/SOURCE_SUMMARY.md)
cover 7,216 cells (44 variants × 164 RNA products), with completed Isovar counts
for 6,202. The other cells explicitly record failed acquisition (442), missing
indexes (220), or unmapped intermediates (352). Five high-depth mitochondrial
cells retain read counts but time out in both modes: nine mode-level timeouts
are during RNA-sequence creation, and one is during independent validation
after Isovar produced 4,154 translations. A missing validated report window
in that last case does not mean Isovar produced no protein.

There is callable alternate evidence for 43 variants and a default protein
window for 39, somewhere among the products. Default outputs contain 749
reference-plus-edit matches and 16 correctly translated RNA differences.
Across completed default-input and primary-only modes, all 83,913 recorded pre-cap ranked windows
pass independent frame/translation/interval/stop checks; these are repeated
window observations across products and settings, **not** independent proteins.

ABCF2, MAP2, MYO15B and MYO9A have alternate evidence but no RNA candidate under
the default coverage requirements. Most have only one allele-read object;
two MYO9A mate alignments can become one object. Counting two alignments is not
proof of two compatible cDNA sequences. ADGRF5 has no callable alternate evidence
in the audited products. NR2F2 has evidence only in the native GRCh37 CeGaT
product, with additional RNA differences. No variant's source-vaccine inclusion
is treated as proof that Isovar ought to reconstruct an expected peptide.

## Scope and source identities

The pinned [variant index](https://osteosarc.com/variants/) contains 44
vaccine-included genomic alleles: 37 SNVs and seven deletions. Exact variant
IDs, not gene names, join membership to the original allele-count export.
There are two DYNC1H1 loci. Nearby MAP2 deletions are not treated as one event.

The original BAM catalogue supplies 50 RNA products. The expanded inventory
reconciles these with the data page, the dated bucket inventory, and paginated
live listings of 17 RNA/analysis-related prefixes. It records all 841 discovered
BAM/CRAM products and each inclusion/exclusion decision. This is not a claim
that imaging/reference-only prefixes were exhaustively re-downloaded.

After inspecting native headers and reconciling DNA catalogue entries, 164
products remain RNA candidates: 156 genomic-coordinate products and eight
unmapped PacBio intermediates. HLA-only, transcript-coordinate, DNA, and
receptor-only alignments are not substitutes for genomic RNA reads. One
uncatalogued `scratch/lms_svaba/all.contigs.bam` remains an unresolved assay:
its genomic header has no program/read-group provenance, so it is not guessed
to be RNA. The complete source-disposition inventory preserves it.

Products are **not independent biological samples**. In particular:

- Tagged/deduplicated ONT products and original/vendor/reprocessed alignments
  can reuse the same underlying material. Published MD5 matches are source
  claims, not verification that we downloaded/hash-checked entire remote BAMs.
- BG009368/SARC0277 have conflicting T0 versus T1/T2 labels in different
  source documents. The matrix retains both claims rather than resolving
  them from expression patterns.
- Blood captures, sample-level outputs, unassigned outputs and enriched
  populations remain separate products. The 67 inspected Cell Ranger configs
  do not establish donor assignments. An observed barcode in a pooled library
  is not automatically this patient's cell, a malignant cell, or a UMI molecule.
- Source variant pages are preserved as claims, not protein oracles. NTF3's
  single-base index entry conflicts with a two-base cDNA change; MT-ND5's
  genomic missense entry conflicts with a `c.774delC` text label.

## Counts, defaults and failures

The full-region matrix is separate from the deliberately selected offline
fixtures. Original indexed queries retain every alignment overlapping the
expanded allele interval, including deletion anchors. The union is fetched
once (`samtools view -M -X`), preserving record multiplicity without duplicate
emission for overlapping queried intervals. No BAM bases, qualities, tags,
CIGARs or flags are rewritten. No full-BAM download or realignment is hidden
inside acquisition.

Every RNA product has a row for every vaccine variant. Unavailable rows have
null counts, not fabricated zeroes. Missing indexes, acquisition failures,
absent genomic coordinates, missing/ambiguous contigs and incomplete audit
work are separate states. A pipeline exception retains its stage and any
successfully completed allele counts. A 120-second audit budget per
source/locus/mode makes `audit_timeout` explicit; it does not change production
support thresholds or fix [#227](https://github.com/openvax/isovar/issues/227).
Timing outcomes depend on hardware/load and are not biological negatives.

The public default returns at most one protein. Audit instrumentation captures
the fully ranked list immediately before that final cap and returns exactly
the configured default slice to the public API. `proteins` is the default
result; `uncapped_ranked_proteins` contains all pre-cap alternatives with
independent checks, including their differing RNA haplotypes. Their separate
validation status must not be mistaken for a default-output failure. No RNA
assembly, threshold, translation or ranking preference changes for this capture.
Deep loci can produce thousands of alternatives sharing the same read names.
The report retains every name once in each row's `supporting_read_name_table`;
each protein's sorted `supporting_read_name_indices` selects its exact subset.
`report.protein_supporting_read_names(row, protein)` decodes and verifies the
subset SHA256. This is lossless storage, not replacement with hashes alone.

The public API's default settings are recorded in each run identity:

- MAPQ at least 1; secondary alignments allowed; duplicates excluded;
  soft-clipped bases excluded; overlapping mates merged.
- Missing base qualities are rejected. There is no new Q20 cutoff or
  probabilistic weighting in this PR; that work remains in #8/#26.
- Balanced context, desired peptide length 25, context ceiling 49 amino acids,
  support fraction 0.85 and minimum RNA-sequence coverage 2; assembly disabled;
  reference prefix at least 10 bases, at most two prefix mismatches.
- Independently validated transcript whitelist, explicitly empty when no
  complete reference model exists. This is a documented isoform restriction,
  not the unrestricted default annotation universe.

Primary-only analysis additionally excludes flags 256, 1024 and 2048 before
running the same API. Exact-CIGAR observations and Q0/Q10/Q20/Q30 sensitivity
are independent descriptive measurements, not calibrated joint confidence.
Q20 reports a nominal 1% base-call error probability, not 99% probability that
the variant is real. MAPQ 255 is unavailable, not an ordinary Phred score.
See the [SAM specification](https://samtools.github.io/hts-specs/SAMv1.pdf).

`reads` counts original supporting alignments via `source_read_count`;
`read_objects` counts post-merge objects; `template_names` counts unique names
within each allele group, and `all_template_names` is their global union.
Compatible protein-supporting names need not span every base in that protein
window. Retained cDNA minimum coverage and mutation-containing peptide-window
counts are reported separately. Distinct `(RG,CB)` and `(RG,CB,UB)` values are
tag observations, not corrected molecules or validated patient assignments.

## Independent protein expectations

Expectations are derived from original Ensembl GTF exons/CDS, cDNA and peptide
FASTA, without using Isovar/Varcode consequence calculations as the oracle.
Every retained complete reference coding model first has to translate exactly
to its original reference protein. Variant reference alleles are checked
against independently retrieved genomic intervals and the original cDNA.
Strand-aware edits then generate transcript-specific mutant expectations.

GRCh38 uses Ensembl 87 (131 complete models; at least one for every cohort
variant). GRCh37 uses Ensembl 75 (140 models; 43 variants); MYO15B has no
validated complete model there. Original Ensembl 75 FASTA identifiers do not
carry `.version` suffixes: release and original archive/subset SHA256 pin
their identity, and missing version fields are **not invented**.

GRCh37 coordinates are independently lifted through the original UCSC chain;
whole-allele, unique ungapped mapping and destination genomic alleles must
validate. A reversed indel needing re-anchoring is rejected, not guessed.
The [chain format](https://genome.ucsc.edu/goldenPath/help/chain.html) uses
half-open coordinates; the source variant identities remain 1-based.

The CeGaT hg19 mitochondrial reference is 16,571 bases, unlike rCRS's 16,569.
Its ND5 allele is at 12,995, not 12,994. A separately identified reference
derivation maps mitochondrial annotation coordinates through the chain and
records two historical genomic background substitutions relative to the
original Ensembl cDNA. It preserves that original biological cDNA/protein;
it does not assert that all hg19 MT sequence equals rCRS.

For each returned window, the oracle checks actual RNA translation, frame,
mutant interval/junction, frameshift, stop and truncation behavior, separately
from equality to reference-plus-nominated-edit sequence. Correct translation
of additional RNA differences can therefore be `different_top_protein`, while
an inconsistent frame/stop/translation is `protein_validation_error`.
Successful local windows do not establish full-length protein reconstruction.
Mutation-containing 25-mer counts explicitly expose insufficient context.

## Mitochondrial interpretation

[#223](https://github.com/openvax/isovar/issues/223) was a real naming bug:
Varcode's `MT` could not find RNA aligned to `chrM`. The fix resolves a unique
alias, prefers an exact name, and rejects ambiguity; name equivalence never
asserts assembly/sequence equivalence.

MT-ND5 is translated with vertebrate mitochondrial
[NCBI table 2](https://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi#SG2):
TGA encodes W, ATA encodes M, and AGA/AGG terminate. The latter correction is
[#224](https://github.com/openvax/isovar/issues/224), supported by direct
[mtRF1 termination experiments](https://pubmed.ncbi.nlm.nih.gov/37141370/).
All 64 codons and initiation-versus-elongation behavior are tested.

The complete bulk MT-ND5 regression region contains 10 ref / 27 alt / 1 other
callable alignments and reconstructs this independently validated 49-aa
A220T window (mutation at index 24):

```text
WDPQQMALLNANPSLTPLLGLLLATAGKSAQLGLHPWLPSAMEGPTPVS
```

Its TGA/ATA codons require table 2; nuclear translation would stop at its
first codon. This demonstrates correct mitochondrial translation **conditional
on mitochondrial origin**, not proof that the reads originate in mitochondria.

NUMTs are **not ruled out**. A published polymorphic numtA includes MT
12361–13227, overlapping this locus
([primary gnomAD study](https://pmc.ncbi.nlm.nih.gov/articles/PMC8896463/)).
All 27 bulk alternate alignments have MAPQ 255 and NH=1, and lie within that
interval. Neither tag proves origin. Long-read diagnostics examine only the
aligned block containing the variant: large intronic `N` spans cannot count
as continuous mitochondrial context. Extending beyond this one published
NUMT does not exclude other/patient-specific NUMTs, chimeras or alignment
error. No competitive whole-genome/mitochondrial/NUMT-decoy realignment was
performed. Matched DNA, validated sample identity, mates/supplementary
alignments and competitive long-read alignment remain necessary follow-up.

RNA allele fraction is neither DNA heteroplasmy nor cancer-cell fraction.
Mitochondrial copy number and within-cell heteroplasmy vary, and RNA adds
expression weighting; see [single-cell primary evidence](https://www.nature.com/articles/s41588-024-01724-8).
Do not apply diploid `2 × VAF / purity` reasoning. The report leaves DNA
heteroplasmy and cancer-cell fraction null. Topiary/vaxrank should retain the
compartment/origin warning and must not infer clonal, tumor-specific vaccine
eligibility from high mitochondrial RNA support alone.

## Offline fixtures and separately labelled stress cases

`corpus/` has 49 original-read cases: all 44 vaccine alleles plus native
GRCh37, hg19-MT, ONT and PacBio mitochondrial examples and native GRCh37
NR2F2 evidence absent from the checked GRCh38 products. Selection favors
allele balance and some protein-supporting names; it is not a random sample
or a full-source VAF estimator. Per-record SAM hashes/ordinals, source region
hashes, original source URLs, BAM/index checksums and original reference
subsets make the fixture provenance auditable. The complete bulk MT region
is retained, not subsampled.

`stress-corpus/` has four separately labelled alleles in four RNA products:
ACSL6/KTN1 insertions, an EPPK1 MNV, and the source-linked NTF3 compound
`AG>GT`. None substitutes for vaccine-index membership. No exact insertion
alternate evidence was found in these four checked sources; this is an
explicit limitation, not positive real-insertion validation.

NTF3's compound allele yields the expected transcript-specific window in
ONT and bulk RNA. Its PacBio region has one ref and one alt exact-CIGAR
observation, both lacking qualities, so current Isovar rejects both; this
is pinned for the #8 missing-quality policy rather than assigning fake scores.

## Reproduction

All acquisition is opt-in and uses bounded network operations. Normal tests
are network-free. Start with a fresh task-specific cache; existing cached
snapshots must match their receipts. A partial/unreceipted artifact requires
inspection or a new cache, never silent overwrite. Do not edit audited source
code while a runner is active: checkpoints pin code, dependencies, settings,
inventory, reference manifest, source BAM and acquisition receipt.

The command-line entry points (each has `--help`) are, in order:

1. `inventory.py CACHE` snapshots membership/catalogue/export metadata.
2. `discover.py CACHE --suffix=-v2` reconciles live listings;
   `metadata.py CACHE` snapshots variant pages, library configs and published
   checksum claims. Header survey uses `acquire.py CACHE --mode headers
   --inventory inventory-live-v2.json`.
3. `liftover.py CACHE` validates original and mapped genomic alleles.
   `references.py ORIGINAL_ENSEMBL_ARCHIVES INVENTORY OUTPUT --assembly ASSEMBLY
   --release RELEASE` builds independently validated reference subsets.
   `cegat_reference.py CACHE GRCH37_REFERENCE OUTPUT` builds the distinct
   hg19-MT annotation derivation.
4. `acquire.py CACHE --mode regions --inventory VALIDATED_INVENTORY
   --assembly ASSEMBLY` fetches original full regions. Nuclear/mitochondrial
   retry partitions use `--retry-failed --subset nuclear|mitochondrial`;
   they remain separately receipted and never reuse failed partial BAMs.
5. `runner.py CACHE REFERENCE OUTPUT --inventory VALIDATED_INVENTORY
   --time-limit 120` audits completed sources. Use separate native GRCh37 and
   CeGaT runs, and `--acquisition-label GRCh38-nuclear` for successful nuclear
   retry partitions. Never mix their coordinate systems.
6. `stress.py CACHE`, then `--observe`; build its references with
   `inventory-stress.json`, and use `--corpus OUTPUT --cache REFERENCE_CACHE`
   to preserve complete original stress regions. `fixtures.py` builds the
   main bounded corpus from explicitly chosen baseline selection rows.
7. `mitochondrial.py CACHE OUTPUT_JSON` generates the limited origin
   diagnostics; `report.py CACHE OUTPUT --audit RUN_DIRECTORY` (repeat the
   audit option for each assembly/partition) assembles the matrix. It refuses
   silently unfinished acquired rows. An explicit diagnostic-only
   `--allow-incomplete` does not certify completion.

Generated reports include checksums. Reassembly from the same checkpoints
is deterministic; fresh timed runs may complete additional protein work on
faster hardware. Compare completed counts/sequences and preserve resource
limits as separate outcomes. Full downloaded RNA stays out of the repository;
only small original-read fixtures and the compressed audit report are kept.

The recorded numerical/protein run source bytes match commit
`122e4a9cb68d1f2f91e9f4f9332af5a499a28eeb` (65 files verified against per-run
SHA256). Reporting/projection changes and optional distribution-metadata
handling are separately pinned in the report-generator manifest. The
`counts-index.json.gz` projection supports routine full-matrix integrity/state
tests without loading millions of support-name indices in each CI worker.
