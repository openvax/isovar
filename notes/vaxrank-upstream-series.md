# Isovar support for Vaxrank

Ship separate reviewed, versioned PRs for the remaining upstream contracts.
Preserve the working checkout and run each release from a clean master.

1. **#589: SV placement support.** Add occurrence-level RNA support computed
   from that occurrence's full-interval, junction-linked witnesses. Preserve
   the candidate union separately and retain the existing v4 export contract
   through additive fields. Test unequal/empty occurrence evidence, overlapping
   read sets, mates, read groups, sample/source scope, lineage, deterministic
   serialization and unchanged input. Exercise the downstream Vaxrank reader.
2. **#571: shared Proteome adapter.** Extend the existing producer adapter for
   SV/fusion and variant-associated ncORF outputs. Preserve source records,
   every alternative/placement, support links and explicit unresolved results;
   never promote an exploratory sequence to observed translation. Use one
   implementation for native and offline inputs with matching API/CLI docs.
3. **#100: patient-background peptide pairs.** Reuse cDNA edit coordinates and
   provenance to construct declared focal-edit or somatic-haplotype baselines,
   retaining supported germline changes. Distinguish RNA-observed baselines
   from computed counterfactuals and unresolved phasing. Cover substitutions,
   indels, frame changes, ambiguous transcripts and incomplete contexts.
4. **#531: original SiHa fixture.** Retain a bounded existing original-RNA
   example with its source/reference hashes and offline reconstruction test,
   then exercise the Vaxrank export contract. This is data retention and
   verification of the existing analysis, not a new discovery algorithm.

File any additional reproducible bugs in the responsible repository and fix
Isovar defects in the appropriate PR. Consumer-specific Vaxrank design/report
bugs (#703/#705) remain separately owned; coordinate a downstream change if an
upstream API change requires it. Do not close broader evidence-acquisition
issues merely because an interface is implemented.

Each PR requires regression tests, lint, the complete test suite, API/CLI/doc
review, green GitHub CI, merge and clean-master `deploy.sh`, followed by live
PyPI/tag verification. Review the next dependency after each release.

Evidence units follow the existing source-scoped RNA identity contract. The
SAM specification distinguishes templates, sequenced segments, supplementary
alignments and alternative placements; alignment rows are not new fragments:
https://samtools.github.io/hts-specs/SAMv1.pdf (sections 1.2 and 1.4).
