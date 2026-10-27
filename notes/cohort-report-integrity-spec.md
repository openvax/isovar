# Cohort report integrity fixes

Fix #596, #597 and #598 in a new PR, then review the complete diff for API/CLI,
documentation and evidence consistency before merge and release.

The review of PR #599 also reproduced #600 (same-event artifacts overwritten
across libraries) and #601 (reference annotation identity collapsed by transcript
ID). Keep fusion artifacts in task-specific directories and resolve reference
comparisons by their recorded identities, preserving ambiguous acceptors.

- Bind rendered results to recorded source, sample, library, task settings and
  input digests. Validate before replacing output manifests or summaries. Permit
  presentation-only edits and relocated input files with the same digests; do not
  require the original BAM during rendering or pretend old results used new code.
  Give additional task families their own receipts, including referenced SV
  artifacts. Reject missing/changed receipts explicitly; never bless existing
  output by hashing it for the first time during rendering.
- Keep donor-anchored fusion proteins renderable when an acceptor coding model
  is unavailable. Project parental codons using the translation start and both
  junction endpoints. Inserted sequence has no parental counterpart; preserve
  partial and out-of-frame interpretations and original coverage.
- Add regression tests using valid reconstruction outputs, source-separated
  manifest rerenders, corrupted artifacts and presentation-only rerenders. Check
  CLI errors and document the same contract exposed by the Python API.
- Render representative fusion PDFs, inspect images, run lint and the complete
  suite, open a version-bumped PR, review it, then require green CI before shipping.

Coordinate semantics checked against the Ensembl core schema (exon phase) and
NCBI's eukaryotic annotation guide (partial CDS/codon_start):
https://useast.ensembl.org/info/docs/api/core/core_schema.html
https://www.ncbi.nlm.nih.gov/genbank/eukaryotic_genome_submission_annotation/
