# One genetic-code source (#682)

All Isovar genetic-code objects should use Varcode's packaged NCBI tables.
Build `standard_genetic_code` through `genetic_code_by_id(1)` and delete the
hand-maintained standard/mitochondrial assignments and implicit wider-start
variant. Keep the `GeneticCode` interface used by SV, CTA, frame and Vaxrank
consumers, including its explicit custom-copy facility.

Retire `translate_cdna`'s unused `mitochondrial` shorthand and require an
explicit keyword-only `translation_table`. Every production caller already
supplies a table except two fusion-report translations; these use the validated
`FusionReference.genetic_code` (currently restricted to table 1). Migrate test
callers and the checked-in fusion/DDX17 report builders, and document the
Python-helper change. No CLI or native evidence schema changes. Vaxrank's
`standard_genetic_code` and `genetic_code_by_id` consumers stay compatible.

The published Vaxrank 3.52.0 source additionally validates annotated ORFs through
`translate_cdna(..., translation_table=None)`, following its resolved #814.
ORFdata's native annotated models and Isovar's own reconstruction use this
nullable contract too. Preserve its standard decoding through the same registry,
without changing undeclared source metadata into an asserted code. A native
TP53 annotated-model regression reproduces the failure when None is rejected.
Check the exact published downstream source, not only an older local checkout.

Non-AUG initiation remains model-specific. Removing an implicit expanded start
set does not remove supplied-ORF initiation handling or change internal codon
decoding. Do not broaden generic fusion/CTA table-1 restrictions. NCBI table 1
has ATG, CTG and TTG starts; alternative model-supported starts remain owned by
their reference-model APIs.

Validation: independent published TCAG-order strings for all 64 codons in
tables 1 and 2, initiation versus elongation, table 3, explicit missing/invalid
table failures, native mitochondrial ORFs, existing source-declared non-AUG
models, and fusion/CTA/SV/viral and Vaxrank consumers. Then lint, full tests,
four-version CI and clean-master release after #689.

Primary references reviewed:
- [NCBI genetic-code tables](https://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi)
- [Lee et al. 2012, global mapping of translation initiation sites](https://pubmed.ncbi.nlm.nih.gov/22927429/)
- Varcode's `effects/codon_tables.py` and packaged `data/genetic_codes.json`.

NCBI's initiator assignments and internal codon assignments are distinct.
The Lee et al. study is not grounds to apply a global expanded start set to
every reference model. This cleanup makes that choice explicit at the caller.
