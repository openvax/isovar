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

For direct translation of a sequence with a known code:

```python
from isovar.genetic_code import translate_cdna

assert translate_cdna("GTGCTGTGAAGATAG", translation_table=3,
                      first_codon_is_start=True) == ("MTWR", True)
```

`mitochondrial=True` remains a compatibility shorthand for table 2. Unsupported
codes and codes requiring context-dependent stop recoding raise an error.
The definitions follow [NCBI's genetic-code tables](https://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi).
Generic exploratory SV ORFs retain their documented standard-code scope; this
does not infer a shared genetic code for mixed-compartment fusion candidates.

## Catalogues

- New code can use `orfdata.load_catalog()` for all 7,264 human GRCh38 ncORF
  models, including lncRNA, processed-transcript, internal and downstream ORFs.
- Other species use explicitly scoped local catalogues; no non-human
  observations are bundled.
- The upstream-only dataset still contains 3,771 models. Isovar retains its
  released `load_uorf_catalog`, `UORFRecord` and initiation-context API as
  compatibility imports, and `isovar uorf-evidence` queries and exports this
  reference. The legacy `annotate_sv_uorf_evidence` call keeps upstream-only
  defaults and field names.

## Osteosarc example

Use `--orf-evidence` in `examples/osteosarc_sv_orfs.py` (`--uorf-evidence` is an
alias). For Sid TPST1–CRCP, 30/61/72 aa alternatives retain 15/1/1 full-interval
RNA witnesses. None retains the native TPST1 MS peptide KIILSKISV. RNA polarity,
initiation, mutant protein expression and mutant MHC presentation remain unresolved.
