# DDX17: retained RNA before frame assessment

The shared `isovar cohort` command now retains local ALT RNA at COLO829 DDX17
`chr22:38489856 GA>G` instead of returning before assembly. In the existing
SRR8615617 selected-site BAM, the replay observes 221 ALT, 325 REF and 56 OTHER
fragments. Default coverage/context settings retain 491 local ALT RNA contexts,
up to 149 nt. These include alternative lengths and placements; they are not
491 distinct molecules or isoforms. The explicit frame state is
`no_coding_reference_context`, and no DDX17 protein is inferred.

GOLGA5 `chr14:92811641 GAGA>G` is a coding deletion control: 109 ALT fragments,
RNA frame candidates retained, and the existing deletion protein/reference
alignment remains available. This two-nomination replay is not a variant census.
The complete source-specific values, settings, module hashes and native result
digest are in [replay.json.gz](replay.json.gz).

From the repository root, with the original pinned BAM available at the path in
[manifest.json](manifest.json):

```sh
isovar cohort --manifest reports/variant_rna/ddx17/manifest.json \
  --output-dir output/ddx17-rna-first
```

The command verifies the BAM and nomination checksums, produces the complete
native JSON/evidence sets, and renders a focused report. Change the manifest's
BAM path when needed, retaining the SHA-256. The source accession is
[SRR8615617](https://www.ncbi.nlm.nih.gov/sra/?term=SRR8615617); the pinned BAM is
the previously aligned, nominated-site subset, not an unrestricted public BAM.
Annotation is Ensembl 87 / GRCh38. No ncORF catalogue or matched-normal evidence
was supplied for this focused assembly regression, so neither is a negative
screening result. No new tumor-specificity conclusion is made.

REF support here is at the allele. This replay does not recover all original
splice-linked reads, reconstruct REF retention paths, resolve cryptic exons or
establish translation. Those #583 requirements remain open. See the
[shared assessment API](../../../docs/variant-read-assessment.md) for the scope
of this step. Visual review of all four report pages also found the existing
plot and absent-catalogue reporting defects tracked in #623 and #624; those
must be fixed before treating this as a finished report.
