# Competing allele interpretations

Different callers and catalogues can describe the same locus differently. One
source reports a single SNV, another a two-base substitution, a third two
adjacent SNVs, and a fourth a nearby deletion. `isovar allele-interpretations`
asks which of these the RNA reads actually carry. Each read's exact sequence
across the locus is compared with each interpretation. Interpretations that give
the same sequence cannot be told apart by any RNA and are reported together.

The result is a table of candidates against evidence. It does not choose an
interpretation or translate one. To get the protein, pass a supported
interpretation's variants to `isovar run` or `run_isovar`. Other variants found
in the same RNA are then reported as co-somatic or germline (see the
[README](../README.md#other-variants-in-the-assembled-rna)).

## Quick start

```json
{
  "reference_name": "GRCh38",
  "interpretations": {
    "catalog A>G": [{"contig": "12", "start": 5494381, "ref": "A", "alt": "G"}],
    "mutect2 AG>GT": [{"contig": "12", "start": 5494381, "ref": "AG", "alt": "GT"}],
    "strelka A>G + G>T": [{"contig": "12", "start": 5494381, "ref": "A", "alt": "G"},
                          {"contig": "12", "start": 5494382, "ref": "G", "alt": "T"}]
  }
}
```

```sh
isovar allele-interpretations --input interpretations.json --bam tumor-rna.bam \
    --sample-id tumor-1 --output interpretations.result.json
```

Each interpretation is a list of variants applied together, with 1-based
`start` as in a VCF. All must be on one contig. The command takes the usual read
options (`--min-mapping-quality`, `--drop-secondary-alignments` and so on), and
`--genome` if the annotation differs from `reference_name`. From Python, use
`reconcile_allele_interpretations(interpretations, alignment_file,
sample_id=..., source=...)` with a dict of interpretation ID to varcode Variants.

## Reading the result

| Field | Meaning |
|---|---|
| `locus.window`, `reference_sequence` | The compared interval (0-based, half-open) and its reference bases |
| `locus.padding` | Exonic bases added on each side (`--flank`, default 10). Padding stays within the exon, so an indel placed elsewhere in a repeat still falls inside the window. It is `[0, 0]` without annotation for the locus |
| `candidates[].haplotype` | The window's sequence under each interpretation |
| `candidates[].rna_status` | See below |
| `candidates[].indistinguishable_from` | Other interpretations with the same haplotype, which no RNA can tell apart |
| `reference_allele` | Support for the reference sequence |
| `observed_alleles` | Every sequence read across the window: its support, the interpretations or reference it `matches`, and its `difference_from_reference` as one trimmed replacement |
| `informative`, `set_aside` | Reads spanning the whole window; segments that overlap it without spanning it, and segments whose alternative placements read different alleles |
| `evidence_sets` | Hashed read IDs for each support, scoped by `[sample_id, source]`, as in the [protein hypothesis export](protein-hypotheses.md#rna-support-and-evidence-sets) |

| `rna_status` | Meaning |
|---|---|
| `supported` | At least `--min-fragments` (default 2) fragments carry exactly this haplotype |
| `contradicted_by_informative_evidence` | No fragment carries it, while at least that many span the window |
| `insufficient_RNA` | Too few fragments either way |
| `unassessed` | The locus could not be compared; `reason` says why, for example unknown reference bases between distant variants |

Only reads that span the whole window count. That makes the comparison exact,
but long windows lose reads. Check `set_aside.segments_not_spanning_window`. An
allele that matches no interpretation is still listed, with its difference from
the reference: the RNA may carry something no caller reported. Absence of RNA
support is not evidence that a DNA variant is absent; the gene may simply not be
expressed.

## Example: the osteosarc patient

Original RNA from the pinned osteosarc fixtures, with calls from the patient's
2024 tumor-versus-normal WGS:

| Locus | Interpretations | RNA |
|---|---|---|
| NTF3 | Catalogued vaccine target A>G; Mutect2 AG>GT; Strelka A>G plus G>T | The two caller forms are indistinguishable and supported by 9 fragments. The isolated A>G is contradicted, and 5 fragments are reference |
| GLIS3 | Catalogued 8-bp deletion; that deletion plus Strelka's A>T; Strelka A>T alone; Mutect2's 7-bp deletion | Only the deletion plus A>T is supported (10 fragments). It is the same haplotype as the single replacement ATGTGGTGA>T at 3856152, which no caller reported |

Both vaccine targets are right about the mutation but miss a second somatic
change on the same RNA. For NTF3 it changes the encoded amino acid (serine, not
the arginine the isolated A>G predicts). For GLIS3 the RNA protein window is
unchanged from the deletion alone.

These are regression tests (`tests/test_allele_interpretations.py`).
