# COLO829 whole-run RNA atlas: implementation plan

Issue: [Isovar #541](https://github.com/openvax/isovar/issues/541).

1. Validate the public SRR8615617 SRA archive and stream its paired reads
   through an aligner against GRCh38, processing the **entire run**.
   Retain indexed original-read alignments at every SMaHT SNV/indel locus and
   both independent SV benchmark breakpoint sets. Record ENA file checksums,
   run size, reference identity, tools, and alignment filtering.
2. Join the original-read evidence back to every source inventory row by
   `(sample, source, class, call_id)`. For all COLO829 SMaHT small variants,
   report reference/alternate/other read and fragment counts, coverage status,
   and translation status. Translate only eligible alternate-supported calls;
   leave unsupported or incompatible cases explicit. Keep TAGET polished
   isoforms in their existing separate tables.
3. For each SV VCF record, report breakpoint coverage and directly observed
   candidate split/pair witnesses separately from any sequence-resolved RNA
   path. Do not infer an ORF from DNA breakpoints or normal splicing. Keep v2
   and lifted-2022 calls distinct, including VCF breakend pairs.
4. Add bounded offline tests for alternate-supported translation,
   reference-only, no-coverage, indel, and an informative SV breakpoint case.
   Compare an independent sample of atlas rows against direct Isovar calls.

The published atlas must state assay/culture mismatch: SMaHT DNA is matched
COLO829/BL, whereas SRR8615617 RNA is a CCLE culture and TAGET RNA is another
source. Whole-run coverage is measured only for this RNA run. ``not_evaluated``
must never be used to imply reference or absence of a mutant isoform.
