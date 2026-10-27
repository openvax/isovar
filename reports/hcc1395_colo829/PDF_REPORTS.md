# HCC1395 and COLO829: variant landscape and RNA-supported candidates

The [HCC1395 report](../../output/pdf/HCC1395_somatic_RNA_evidence.pdf) and
[COLO829 report](../../output/pdf/COLO829_somatic_RNA_evidence.pdf) include
inventory and support landscapes, category exemplars, and focused protein
alignments. The PDFs have bookmarks for the cohort evidence and deep dive.

Standalone deep dives:

- [COLO829 CDKN2A](focused_proteins/COLO829_CDKN2A_frameshift.pdf): 5,248
  fragments for the reported 24-aa window, the full reference protein for
  orientation, highlighted mutant tail and an assembly/reference track.
- [HCC1395 EIF3K-CYP39A1](focused_proteins/HCC1395_EIF3K_CYP39A1_fusion.pdf):
  531-aa coding candidate, both parent references, junction zoom, exact
  fragment support per codon, and coding-path assembly coverage.
- [Candidate protein FASTA](focused_proteins/candidate_proteins.fasta) keeps
  the partial CDKN2A window and complete fusion candidate explicitly separate.

Every library, DNA nomination set, polished input and reference is listed in
[SEQUENCING_SOURCES.md](SEQUENCING_SOURCES.md). Source-specific counts must not
be pooled. These reports cover the selected inventories, not every mutation
in every culture of either cell line.

## Organization

1. DNA-to-RNA small-variant funnel and **parallel SV discovery routes**.
   Split alignments, cross-end mates, and soft-clipped/unmapped sequence are
   complementary. The old split/mate screen is not an eligibility filter;
   its unenumerated soft-clip channel is unknown, not zero. Exon overlap is
   an annotation, not a requirement for reconstruction.
2. Inventory composition and a log-scale landscape of ALT-fragment support
   against support for the top local coding window. Baseline points lack a
   reconstructed local context under the recorded settings.
3. Full source-record accounting, RNA-support distributions, the explicitly
   exon-restricted legacy SV screen, and original-read junction replay.
4. Up to **three ALT-supported examples per category**: high-support SNVs,
   frameshifts, other indels, SV-derived RNA transcripts, and ncORF overlaps.
   All examples include gene names and copyable one-based IGV coordinates.
   Fewer than three qualifying examples remain an explicit result.
5. Reference-versus-RNA comparisons, remaining category gaps, sources, and
   the appended focused CDKN2A or EIF3K-CYP39A1 protein report.

Only differing focal mutant residues use the focal-change color. Reference
residues stay neutral. Germline and other phased somatic changes require
source-qualified evidence; unexplained RNA differences are marked unclassified.
For example, the HCC1395 TPI1 DNA call is synonymous but its reported RNA
window differs outside the focal site. Its comparison is anchored by the
reference genomic codon and does not establish the origin of those differences.
A reference-identical or reference-incompatible context is not silently counted
as a validated neoORF; broader classification remains in [#551](https://github.com/openvax/isovar/issues/551).

In-frame deletions now have a verified reference gap and a vertical focal-color
marker; unchanged residues on both sides remain aligned and neutral. GOLGA5 is
on COLO829 page 12. Its original `p.K404del` label has an upstream
[Varcode numbering defect](https://github.com/openvax/varcode/issues/643); the
independently verified one-based deleted residue is K405.

COLO829 page 13 adds the [DDX17 original-RNA audit](ddx17_audit/README.md):
27,309 recovered pairs, exon/intron coverage, alternative placements and explicit
reference-frame hypotheses. The original coding-reference gate did not assess
cryptic splicing. Intron-containing RNA and alternative junctions are present,
but no mutation-phased coding path is established.

Protein-window support, focal allele counts, junction support, local branch
support and whole-path tiling are different measurements. The complete HCC
fusion is assembled from RNA-witnessed bases; the CDKN2A report does not fill
its missing flanks from reference. No RNA reconstruction measures translation
or antigen presentation. Mutated ncORF proteins and viral/other antigen-source
surveys remain unassessed or unresolved where stated.

## Historical benchmark denominators

The exon-restricted screen below is retained for provenance. It is not the
new SV discovery funnel and does not contain the separately supplied HCC fusion.
The [unified COLO829 audit](unified_sv_audit/REPORT.md) accounts for all 188 SV
source records and all retained paths without an exon-overlap eligibility gate.

| Measure | HCC1395 | COLO829 |
|---|---:|---:|
| Ingested small variants | 41,072 | 46,064 |
| Ingested SV source records | 1,788 | 188 |
| SV source records involving an Ensembl 87 exon | 354 | 97 |
| Exonic SVs with RNA alignment near either nominated end | 172 | 64 |
| Exonic SVs with split-read or cross-end mate candidates | 6 | 6 |
| Exonic SVs with an original-read RNA junction in selected replay | 3 | 0 |
| Any Ensembl 87 exon involved | 1,728 | 1,395 |
| Exonic with original RNA ALT | 292 | 502 |
| Exonic with a local protein context | 22 | 124 |
| Ensembl 87 coding-sequence effect | 491 | 373 |
| Pinned orfdata ncORF coding-block overlap | 17 | 5 |
| Either coding annotation | 505 | 378 |
| Original RNA ALT-observed calls | 343 | 1,933 |
| ALT-observed calls with a local protein context | 22 | 124 |
| Ranked local protein contexts | 119 | 12,218 |
| DNA-predicted top-effect frameshifts | 9 | 5 |
| Resolved SV-derived proteins | 0 | 0 |


## Reproduce and inspect

The pinned reference/locus tables keep rendering offline. To refresh them from
the local Ensembl 87 cache, use `build_report_locus_annotations.py` and
`build_example_reference_contexts.py`. Then render:

```sh
python -m reports.hcc1395_colo829.build_focused_protein_reports
python -m reports.hcc1395_colo829.build_evidence_pdfs
```

The [locus annotation table](report_locus_annotations.tsv) includes all
ALT-supported small variants and all benchmark SV source rows. Protein genes
come from compatible reference transcripts; other genomic overlaps and SV-end
genes remain separately identified. Unresolved partners are never assigned an
invented IGV location. The [reference-window table](example_reference_contexts.tsv)
pins sequence, transcript, residue offset, reference-protein digest and placement
method. Focused report JSON files additionally pin every input and the exact
per-codon support used for the figures.
