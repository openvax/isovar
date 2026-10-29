# RNA-first literal-variant assessment (#583, #622)

The current translation list API returns before assembly when no coding context
exists. Preserve its return type, but delegate to a typed assessment that always
retains allele-supporting RNA or a specific assembly failure. Keep coding-path
compatibility safeguards: compatible groups and unmatched RNA use the same
assembler; unmatched reads must not be forced into an annotated CDS. Avoid
mutable last-result state on a creator. Existing ranking should consume the
assessment rather than recompute it, and IsovarResult/export/cohort should expose
RNA candidates and the explicit reason a frame was not assigned.

Use the same flow for coding, noncoding and intronic literals. Reuse the supplied
ORF interface for reference models and keep ordinary transcripts as competitors.
Fix the remaining _context string replacement (#622) with one reused Varcode
resolver over the supplied reference vocabulary. Exact custom/alt/viral/decoy
names win; preserve missing and ambiguous resolution rather than equating names
or counting an unresolved chromosome as unsupported RNA.

Stage the work into reviewable PRs if needed:
1. Typed assembly/translation assessment; no early annotation exclusion; shared
   CLI/API/cohort exports, no-coding-context reason; #622; existing coding outputs
   and Vaxrank compatibility. Reproduce original DDX17 local RNA and REF controls.
2. Evidence-driven original-read recovery and placement-aware RNA paths. Reuse
   original-library acquisition, read-end policy, allele identities and shared
   exon/CDS anchor helpers from the SV path workflow. Keep each library separate.
   Compare retained-intron/cryptic-exon/ordinary-reference alternatives, annotate
   allele-to-junction linkage, support partial or complete frame candidates only
   when the observed path connects an upstream coding anchor or explicit ORF.
   ATG-only hypotheses stay exploratory. No reference-filled sequence is observed.
3. Run the same CLI on DDX17 and cohort inputs, report ascertainment/assembly
   budgets, normal evidence and unresolved states. Keep #583 open until original
   recovery, linkage and frame/no-frame outcomes are reproducible through the
   library workflow, not just an extra report-script scan.

Primary sources checked 2026-10-29:
- Smart et al., 2018, https://pmc.ncbi.nlm.nih.gov/articles/PMC6226333/ and
  https://www.nature.com/articles/nbt.4239: retained introns can produce
  experimentally observed peptides; depth or retention alone does not establish
  translation of a particular sample's RNA or mutation-induced splicing.
- https://useast.ensembl.org/info/genome/genebuild/biotypes.html and
  https://useast.ensembl.org/info/genome/genebuild/transcript_quality_tags.html:
  retained-intron/noncoding and incomplete CDS annotations are distinct from a
  supplied complete coding frame. Do not convert a missing annotation into a
  biological absence claim or invent a start codon.

Current original DDX17 target: COLO829 GRCh38 chr22:38489856 GA>G, deleted A at
38489857, 223 ALT observations / 221 fragments in the prior selected source.
Original audit folder is /Users/iskander/code/isovar/output/ddx17_ref_alt_followup_20261027/original_audit/.
Those are source-specific prior counts, not a newly verified full-library result.
