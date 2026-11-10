# Reference ORFs and translation evidence

Isovar can annotate the exploratory ORFs from [`isovar sv-rna`](sv-rna.md) with
matching reference ORF models and their translation evidence. The reference
catalogues live in [orfdata](https://github.com/openvax/orfdata), alongside
PyEnsembl. See the orfdata README and manifests for primary sources and
limitations.

## Annotate SV ORF candidates

`annotate_sv_orf_evidence(result, transcripts, catalog=...)` binds each matched
reference ORF to the supplied PyEnsembl transcript and delegates modifications
to Varcode's `annotate_orf_product`. The annotation preserves separate reference
model/evidence and the conditional altered product, with candidate-specific RNA
support. JSON and TSV exports preserve this structure and catalogue provenance.
The generic API uses the field names `reference_orf_annotations` and
`orf_catalogue_metadata`.

## Explicit genetic codes

Bound ORFs retain their declared NCBI `translation_table` through Varcode's
variant/RNA-product annotation and Isovar's `annotate_sv_orf_evidence` overlay.
The pinned OpenProt IP_306399 example uses vertebrate mitochondrial table 2;
the native yeast OLI1 and Q0092 fixtures use yeast mitochondrial table 3.
Reference evidence remains attached to the reference product. A reconstructed
RNA sequence does not establish mutant translation or initiation.

Small-variant RNA reconstruction also carries `ReferenceContext.translation_table`.
Contexts with different codes remain distinct even when their local sequences
and frames match. Selection uses Varcode's shared transcript/source declarations
and reference species/compartment metadata; an unknown mitochondrial reference
requires an explicit assignment. The code is recorded in each translation's
`reference_context` in the protein-hypothesis export and used to validate it.

Direct sequence translation defaults to the standard code (NCBI table 1),
including mammalian nuclear translation. Supply another code when needed:

```python
from isovar.genetic_code import translate_cdna

assert translate_cdna("ATATGAAGG") == ("I", True)
assert translate_cdna("ATATGAAGG", translation_table=2) == ("MW", True)
assert translate_cdna("GTGCTGTGAAGATAG", translation_table=3,
                      first_codon_is_start=True) == ("MTWR", True)
```

`translate_cdna` accepts an optional `translation_table` keyword. Use
`translation_table=2` in place of the removed `mitochondrial=True` shorthand;
use table 1 for standard-code translation. Explicit `None` retains the existing
standard decoding convention for annotated ORFs with undeclared source codes,
without rewriting their metadata to assert table 1. This also preserves the
released Vaxrank validation contract ([vaxrank #814](https://github.com/openvax/vaxrank/issues/814)).
Unsupported explicit IDs and codes requiring context-dependent stop recoding
raise an error. This sequence-only helper does not infer species or cellular
compartment; transcript-aware reconstruction resolves those from reference
context. Native ORF reconstruction uses that resolved code for the candidate,
reference window and focal-allele counterfactual, even when the source's
declared code is None.
The definitions follow [NCBI's genetic-code tables](https://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi).
The public `standard_genetic_code` object uses the same Varcode registry as
explicit-code translation. Isovar keeps no separate codon assignments or
implicit expanded-start table. `first_codon_is_start=True` uses the selected
table's initiator set at the first position only. Source-supported non-AUG
starts outside that set remain conditional interpretations of their supplied
ORF model, rather than a change to the global genetic code.
Generic exploratory SV ORFs retain their documented standard-code scope; this
does not infer a shared genetic code for mixed-compartment fusion candidates.

Native `reconstruct_variant_orfs` exports use `isovar.variant_rna_orfs.v2`.
Every protein candidate carries an integer `translation_table` and
`translation_table_source` (`orf_model` for an explicit source declaration,
otherwise `reference_context`). The model's `reference.translation_table`
remains unchanged. The effective code also participates in candidate identity.
Readers should use the shared public contract:

```python
from isovar import variant_orf_translation_table

table = variant_orf_translation_table(candidate, model["reference"],
                                      schema=result["schema"])
```

V2 rejects absent/unsupported effective codes and conflicts with explicit source
declarations or provenance. V1 replay retains the original nullable source-code
convention and does not reinterpret old mitochondrial results. Isovar's native
proteome reader and Vaxrank 3.54+ consume the same contract.

A NuMT is nuclear sequence derived from mitochondrial DNA. If a supplied nuclear
transcript contains such sequence, its nuclear code applies (table 1 in humans);
similarity to mitochondrial sequence or a mitochondrial-looking gene name does
not select table 2. A code assignment is conditional on the supplied reference
placement and does not prove read origin or that the locus encodes a protein.
Short-read cross-mapping can obscure origin; longer context and paired reads
can improve placement ([Maude et al. 2019](https://doi.org/10.3389/fcell.2019.00201)).
Non-reference NuMTs can also create false mitochondrial calls
([Laricchia et al. 2022](https://pmc.ncbi.nlm.nih.gov/articles/PMC8896463/)).
Isovar's osteosarc audit has limited original-read diagnostics and explicitly
reports NuMT origin as not ruled out; it does not perform a general competitive
nuclear/mitochondrial realignment or NuMT screen. A confident NuMT assignment
requires supporting placement evidence, such as reads spanning unique nuclear
flanks, with ambiguous alternatives retained rather than assigned by code choice.

## Catalogues

- New code can use `orfdata.load_catalog()` for all 7,264 human GRCh38 ncORF
  models, including lncRNA, processed-transcript, internal and downstream ORFs.
- Other species use explicitly scoped local catalogues; no non-human
  observations are bundled.
- The upstream-only dataset still contains 3,771 models. Isovar retains its
  released `load_uorf_catalog`, `UORFRecord` and initiation-context API as
  compatibility imports, and `isovar uorf-evidence` queries and exports this
  reference. It prints a table by default; `--json` writes the complete JSON
  document and `--format fasta`, `tsv` or `gtf` export every matching record.
  An unknown ORF ID or a REF that disagrees with the genome exits 1; a rejected
  command line exits 2. The legacy `annotate_sv_uorf_evidence` call keeps
  upstream-only defaults and field names.

## Osteosarc example

Use `--orf-evidence` in `examples/osteosarc_sv_orfs.py` (`--uorf-evidence` is an
alias). For Sid TPST1–CRCP, 30/61/72 aa alternatives retain 15/1/1 full-interval
RNA witnesses. None retains the native TPST1 MS peptide KIILSKISV. RNA polarity,
initiation, mutant protein expression and mutant MHC presentation remain unresolved.
