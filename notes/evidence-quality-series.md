# Evidence quality and unified reconstruction series

Implement the requested issues as reviewed, independently versioned releases.
Keep native API, CLI, JSON/tabular exports and reports on shared library paths.
Preserve working checkouts, original inputs and completed scientific receipts.

1. **Contig resolution:** adopt released Varcode `ContigResolver`, reuse each
   BAM/reference vocabulary, prefer source-exact names and preserve accession
   versions. Expose missing/ambiguous/unavailable contigs separately from measured
   zero support. Test canonical, alternative, viral, decoy and arbitrary names,
   mixed reference contexts, source spelling, ambiguity and repeated access.
2. **#590 / #551:** centralize coding-window support and reference comparison;
   distinguish assembly contributors, exact paired-template coverage and original
   segment spans. Count verified amino-acid changes independently of annotated
   mutation intervals and retain unresolved comparisons.
3. **#595:** inventory matched-normal sources by stable public identities and
   processing relationships. Keep called ALT/reference, no-call, variant-only
   absence and unavailable evidence distinct. Carry this scope into backgrounds
   and cohort provenance without silently revising completed receipts.
4. **#593:** prepare and reuse an original-template index per source/policy.
   Query complete templates without changing discovery/recount eligibility or
   evidence units; demonstrate equivalence and actual inventory performance.
5. **#302 / #583:** audit existing read-end annotation first. Complete shared
   technical-sequence handling and an explicit intronic RNA assessment path,
   retaining biological/technical ambiguity and linking allele observations to
   splice/retained-intron models before proposing coding consequences.
6. **#575 / #588:** verify original source conventions and validation metadata,
   preserve original records, and import only justified oriented adjacencies.
   Unsupported definitions remain explicit unresolved records; never fabricate
   coordinates or treat unassessed calls as zero support.

For each stage: reproduce the problem, write a bounded spec, implement with
clear names/type signatures/docstrings and centralized helpers where useful,
review the final diff as a maintainer, add meaningful regression tests, run
lint/full tests, pass GitHub CI, merge, run clean-master deploy.sh and verify
published artifacts. File newly discovered defects on the owning repository.
Close issues only when their remaining acceptance criteria are fulfilled.

Scientific reference semantics follow SAM v1.6: @SQ SN is the actual alignment
reference vocabulary; alternate-locus AH does not establish identical sequence
or coordinates. Name conversion is not liftover or sequence identity:
https://samtools.github.io/hts-specs/SAMv1.pdf. The released Varcode resolver owns
spelling conversion and ambiguity; Isovar owns input-context lifetime and RNA
assessment/reporting state.
