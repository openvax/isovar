# CTA and viral workflow unification

The current viral and antigen-isoform APIs instantiate fixed read collectors.
Cohort manifests can record a `read_collection` policy that these tasks ignore.
This violates the recorded evidence contract and prevents shared library reuse.

1. Route literal RNA paths, antigen isoforms and viral paths through one default
   read-policy factory. Expose the existing ReadCollector, original-library index
   and custom eligibility identity in both antigen APIs. Preserve defaults and
   retain the effective policy in results. Dedicated CLIs use common RNA flags;
   cohort tasks reuse one validated original-library context per source/policy.
   Test discovery/support agreement, technical sequence settings, cache identity,
   CLI/API/cohort parity and unchanged default evidence. Ship as the #640 PR.
2. Add versioned OncoRef CTA target selection to the existing isoform workflow.
   Keep catalogue selection separate from transcript reconstruction and vaccine
   admission. Pin the provider version, policy, source snapshot and all gene
   mappings; retain missing/ambiguous annotation. Do not copy a gene list or
   maintain a second CTA definition in Isovar. Keep explicit genes available.
3. Inventory available original RNA libraries for Sid, COLO829 and HCC1395.
   Run the shared cohort CLI with catalogue-selected CTA targets and the viral
   contigs actually present in matching references. Preserve library/culture and
   input scope. Where references lack viral placements, recover original RNA
   through explicit competitive alignment before interpreting support. Record
   unavailable inputs and incomplete screens. Generate native results and PDFs
   before the final broader verification, then review figures and sequences.

Existing user checkouts and the shared environment are preserved. Use isolated
feature worktrees and release environments; each PR requires lint, full tests,
green CI, merge, deployment and verification of published artifacts.
