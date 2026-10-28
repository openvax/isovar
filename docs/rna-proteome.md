# Proteomes from RNA hypotheses

`isovar.from_rna_result(result, reference=...)` adapts a native JSON-ready result
or its saved JSON into an ORFdata `Proteome`. It retains every supplied protein
candidate, original source metadata and unresolved outcome. It does not rerun
assembly, fill sequence from reference, select a winning frame, or add evidence
of protein expression.

| Input schema | Mapping |
|---|---|
| `isovar.protein_hypotheses.v2` | Existing small-variant hypothesis adapter |
| `isovar.sv_rna_orfs.v4` | Every exploratory nucleotide candidate and occurrence |
| `isovar.sv_rna_candidates.v5` | Annotated-frame translations, all retained paths and the shared exploratory ORF mapping |
| `isovar.fusion_rna.v3` | Every supplied-path translation, including partial/ambiguous models |
| `isovar.variant_rna_orfs.v1` | Every retained model/RNA/protein candidate, including ncORFs |

The v5 SV route preserves annotated-frame proteins and paths without an ORF in
addition to the compact exploratory export. Its source record links to the
exploratory export record; these are representations of the same RNA evidence,
not independent support. Historical v4 files cannot supply native path data or
annotated-frame products that were never exported. Save the full v5 result to
replay the complete native mapping offline.

```python
import json
from orfdata import ReferenceContext
from isovar import from_rna_result

reference = ReferenceContext(9606, "human", "GRCh38", "Ensembl", "87")
with open("sample.orfs.json") as handle:
    rna = from_rna_result(json.load(handle), reference=reference)
rna.save("sample.proteome.json.gz")
rna.write_fasta("sample.proteins.fasta")

# Merge a separately sourced reference proteome or another RNA source:
# augmented = reference_proteome.merge(rna)
```

For native ordinary small-variant `IsovarResult` objects, the existing
`to_proteome(results, sample_id=..., source=..., reference=...)` API remains
available. Neither route downloads annotation or requires the original BAM
when converting saved results.

The equivalent CLI accepts JSON or gzipped JSON. Save `reference.to_dict()` to
`reference.json`, or provide the corresponding ReferenceContext JSON fields:

```sh
isovar proteome --input sample.orfs.json --reference-context reference.json \
  --output sample.proteome.json.gz --fasta sample.proteins.fasta

isovar proteome --input library-A.orfs.json --input library-B.orfs.json \
  --reference-context reference.json --output combined.proteome.json.gz
```

The standalone `isovar-proteome` command has the same options. Inputs are
validated before output is written; outputs cannot overwrite an input or the
reference-context file. No additional installation extra is needed.

## Evidence and identity

The snapshot retains event, hypothesis, model, RNA, occurrence, translation and
evidence records as appropriate to the source schema. Links retain array order
and source relationships. Original candidate/model/sequence IDs remain in their
native payloads. New RNA-adapter record IDs are qualified by a canonical export
digest so different assessments of the same source cannot overwrite one another.
Producer sequences also carry that digest and the original sample/source IDs.
Repeated identical conversions merge idempotently. Equal amino-acid sequences
from distinct candidates or libraries remain separate by default.

SV occurrence support belongs to that placement. The candidate union is a
separate record and is never copied onto an unsupported or unmeasured placement.
Historical v4 occurrences without scoped support retain that absence and receive
no invented support link. Read sets may overlap; support counts are not additive.
Cell/UMI completeness, signal ancestry, reference comparisons, frame/start
ambiguity, search limitations and witness coordinates remain in native payloads.

Fusion RNA retains original observation records and their library/source IDs;
its producer run ID comes from the supplied assembly provenance's `source`.
Variant ORFs retain reference-model evidence separately from mutant RNA evidence,
along with coding coverage and the focal-allele counterfactual. No RNA-derived
entry gains ORFdata protein-expression evidence. Start/stop completeness describes
a sequence hypothesis, not observed initiation or translation.

Recorded builds must agree with the supplied genomic ReferenceContext. Complete
variant-ORF reference contexts must also agree, including annotation and reference
pins. Where older schemas record only an assembly name, species and annotation
remain explicit caller assertions, not independently verified facts. Missing
producer versions remain unknown; the converter does not substitute its own
installed version for historical reconstruction software.

RNA evidence scopes, set identities, counts and links are checked for consistency.
These checks detect contradictory records; they do not authenticate an externally
supplied result. Partial, zero-product and limited-search outcomes remain visible
even when the corresponding FASTA has no entry.
