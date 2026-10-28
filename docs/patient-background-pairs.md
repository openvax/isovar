# Patient-background peptide pairs

`export_patient_background_pairs` returns every retained local protein window
paired with an explicit somatic-reversion baseline. It operates on native
`run_isovar` results and uses transcript/cDNA edit coordinates to retranslate
the observed background. Supplied germline substitutions and indels remain.
No amino acid is replaced directly, and no missing flanking sequence is filled
from annotation.

```python
from isovar import run_isovar, export_patient_background_pairs

results = run_isovar(somatic_variants, alignment_file,
                     germline_variants=normal_variants)
pairs = export_patient_background_pairs(
    results, sample_id="patient-1", source="tumor-RNA", baseline="focal")
```

The equivalent CLI uses the same pipeline and export. All reconstructed protein
hypotheses are retained by default; an explicit protein-count cap is recorded.

```sh
isovar patient-background-pairs --vcf somatic.vcf --germline-vcf normal.vcf \
  --bam tumor-rna.bam --sample-id patient-1 --source tumor-RNA \
  --baseline focal --output pairs.json.gz
```

The standalone `isovar-patient-background-pairs` command has the same options.
The JSON/JSON.gz output uses `isovar.patient_background_pairs.v1`; it requires
no extra dependencies, annotation downloads beyond the ordinary reconstruction,
or new BAM pass. The original reference-allele observations already collected by
Isovar supply the separate baseline support assessment.

## Choosing a baseline

| Policy | Removed | Retained |
|---|---|---|
| `focal` (default) | Nominated focal somatic edit | Germline, co-somatic and unexplained RNA edits |
| `supplied_somatic_haplotype` | All edits attributed to supplied somatic calls on that RNA | Germline and unexplained RNA edits |

Each record lists the actual removed and retained edits, including transcript
coordinates, their positions in the observed RNA, restored reference bases and
source event IDs. A supplied-variant registry retains genomic coordinates and
declared origins. Equal proteins from different transcripts/nucleotide contexts
remain separate, linked to the standard protein-hypothesis and translation IDs.

These are **bounded coding windows**, not necessarily complete proteins or fixed
MHC peptide lengths. The baseline uses the mutant's recorded RNA coding extent,
adjusted for reverted indels, and ends at its own first in-frame stop. Mutant and
baseline lengths may differ. Incomplete trailing codons and absent stops remain
explicit. A stop removed by reversion cannot justify a downstream sequence that
falls beyond that recorded extent. `reference` stores the local annotation cDNA
and transcript interval used to validate edit replay; it is not patient RNA.

## Evidence and unavailable pairs

`construction: transcript_cdna_somatic_reversion` always identifies a computed
counterfactual. `observation: reference_allele_rna_span` is added only when separate
original reference-allele observations exactly span its translated coding window;
otherwise it is `computed_only`. Coverage, allele-observation spans and single
original-read spans are reported separately. Reference reads incompatible with
the selected transcript are excluded. Neither observation state establishes
normal-tissue expression, the entire returned cDNA, translation, or immunogenicity.

Mutant phase is `single_original_read_span`, `allele_observation_span`,
`assembly_only`, or `incomplete_rna_coverage`. Tiling does not demonstrate a single
molecule. The computed sequence remains conditional on its RNA assembly and frame;
uncertain phase is not converted into a normal genotype claim. Original hashed
read/fragment identities have the same sample/source scope as other Isovar exports.
Support sets may overlap and their counts must not be added.

Exact edit replay is required. Overlapping/out-of-bounds edits, incomplete reference
matching, a reverted edit crossing a coding boundary, or an unavailable focal
coding interval produce `state: unavailable`, a reason and no baseline sequence.
An empty baseline and an event without any protein hypotheses are also retained
explicitly. Legacy reads with unavailable identities retain unknown evidence IDs.

Vaxrank can use these records to offer an explicit patient-background comparator.
Its annotation-reference comparison remains a separate policy until the consumer
selects and checks the paired record's construction, phase and evidence state.
