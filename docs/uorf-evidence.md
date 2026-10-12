# Human uORFs with translation evidence

Isovar packages an offline, versioned reference catalogue of **3,771 human
upstream ORFs**, including upstream ORFs that overlap the canonical CDS
(`uoORF`). Each entry has its reported transcript, GRCh38 spliced coding blocks,
exact reference protein sequence, ribosome-profiling studies, and available MS
peptide observations. It is an evidence reference for investigating variants in
5′ UTRs, rather than a list of validated mutant antigens.

The dataset (`2026-10-08.2`) has 1,219 models with supporting reported MS
observations, including 1,059 with HLA peptides and 31 with clearly classified
shotgun reports. Categories overlap. Legacy mixed-assay reports remain
`ms_unclassified`. There are 16 models with compatible peptides receiving
convincing community spectrum reviews. Those reviews can support a peptide
sequence without identifying which compatible ORF produced it.

## Python API

```python
from isovar import load_uorf_catalog

catalog = load_uorf_catalog()  # Packaged data; no network access.
tpst1 = catalog.get("c7riboseqorf56")
print(tpst1.protein_sequence)
print(tpst1.peptides)

# Genomic positions and coding blocks are 1-based inclusive.
overlaps = catalog.query(
    contig="chr7", start=66240380, end=66240380, genome_build="GRCh38")
annotation = tpst1.map_genomic_position(66240380, genome_build="GRCh38")
print(annotation["amino_acid_position"])  # 32, in the observed native peptide.
print(annotation["reference_ms_support"])  # True; reference observation only.
print(annotation["mutant_translation_evidence"])  # "not_assessed"

# Require convincing reviewed spectra; source uniqueness is a separate filter.
reviewed = catalog.query(evidence="hla", require_reviewed=True)
source_unique = catalog.query(evidence="hla", require_unique=True)
```

`UORFCatalog.query` combines exact gene-symbol/Ensembl-gene, transcript, biotype,
assay and genomic filters. Transcript IDs are the source's **unversioned** IDs;
versioned or alternative IDs do not silently inherit evidence. Genomic queries
require an assembly and intersect actual coding blocks, so an intronic position
does not count as an ORF overlap. Minus-strand and exon-spanning codons map in
transcript order. The final three bases are the model's stop codon; protein
sequences omit `*`. The initial codon identifies the reported model's start,
not an independently confirmed initiation site.

The immutable `UORFRecord` and `UORFPeptide` records serialize with `to_dict()`.
`catalog.metadata` returns a copy containing the dataset version, study PMIDs,
source hashes, licenses and evidence scopes. The packaged file is verified
against its manifest on load. To pin a local compatible JSON/JSON.gz catalogue:

```python
catalog = load_uorf_catalog("reference.json.gz", expected_sha256="<your-pin>")
```

## Initiation context and mutations

Every current model has a reference initiation context from checksum-pinned
Ensembl 101 transcripts. The entire spliced ORF path and protein are validated
before flanking sequence is attached. Context provenance retains the resolved
transcript version. No adjacent intronic sequence or alternative isoform is
silently substituted. Complete/partial/unavailable status and failure reasons
are explicit. This reference release supplies historical model contexts,
rather than proving the patient's transcript or its biological 5-prime end.

```python
from isovar import annotate_initiation_context, load_uorf_catalog

# Also works on reconstructed mutant RNA or partial read/assembly sequence.
context = annotate_initiation_context("ATGG", 0)
print(context.assessment)  # Partial; +4 G observed; cap distance unknown.

native = load_uorf_catalog().get("c7riboseqorf56")
print(native.initiation_context.assessment["display"])
# GAGGCCAGG[ATG]CCGTCC
print(native.annotate_initiation_snv(
    66205480, "A", "C", genome_build="GRCh38")["sequence_changes"])
# ['kozak_preference_lost']; -3, no coding-residue change.
```

`InitiationContext` is immutable and serializes raw sequence, spliced genomic
placements, source identity and its assessment. Preferred positions are
individually true/false/null; missing and ambiguous bases remain unknown.
The -3/+4 counts provide a loose interpretation, with no exclusion threshold
or calibrated initiation probability. Non-ATG and zero-offset starts are
retained. `five_prime_complete=True` requires caller evidence of the biological
5-prime end; only then is `cap_distance_nt` populated.

`map_initiation_position` distinguishes context/start/coding roles. The SNV
API takes single forward-genomic REF/ALT bases, validates REF and converts
minus-strand alleles. A +4 SNV can have both context and coding effects.
Its alternate context and codon consequence are predictions, not mutant
translation evidence. `query(include_initiation_context=True, ...)` includes
placed spliced flanks; the default continues to query coding blocks only.
Neither includes intervening introns. Old v1 catalogues without contexts
remain loadable, with unknown context rather than inferred sequence.

See [the mutation annotation guide](uorf-mutation-annotations.md) for the
identity/location/consequence/evidence organization and the Varcode proposal.

## Command line and dataset exports

```sh
isovar uorf-evidence --gene TPST1
isovar-uorf-evidence --evidence shotgun --format tsv > shotgun-uorfs.tsv
isovar-uorf-evidence --evidence hla --require-reviewed --format fasta > reviewed.fasta
isovar-uorf-evidence --position chr7:66240380 --genome-build GRCh38
isovar-uorf-evidence --position chr7:66205480 --genome-build GRCh38 --ref A --alt C
```

JSON retains provenance and nested evidence. TSV uses explicit JSON columns for
coding blocks, study IDs, peptides and initiation contexts, so exports retain those relationships.
FASTA exports exact reference proteins. All commands work offline after an
ordinary Isovar install. The dataset and its source manifest are installed under
`isovar/data/uorf-evidence/`; the builder and this guide are also in the source
distribution.

## Interpreting evidence and mutations

| Observation or event | What to do with it |
|---|---|
| Ribosome profiling | Prioritize a reported translation model; inspect initiation, reading frame and the patient's expressed isoform. It does not establish stable protein abundance or HLA presentation. |
| Shotgun MS peptide | There is a reported native protein fragment. Inspect peptide mapping and review quality; it does not establish the patient's mutant sequence. |
| HLA MS peptide | There is a reported presented native peptide sequence. Record the source assay and HLA context separately from predicted mutant binding. |
| SNV or indel inside a uORF | Validate genomic REF, reconstruct the patient's RNA isoform and translate in this uORF's frame. Check whether the changed residues overlap the observed native peptide. |
| Start/stop change | Reconstruct initiation loss, extension or truncation explicitly. Reference initiation evidence cannot prove use of a new start or the extended mutant tail. |
| Splice change or fusion | Reconstruct the new splice path and coding frame. Keep evidence for retained native residues separate from the new junction or tail. |
| UTR change outside a catalogue coding block | Investigate regulation, alternative initiation and other isoforms; the catalogue does not establish a protein change. |
| No entry or no MS observation | Leave translation unassessed. Limited study coverage is not a negative result. |

`map_genomic_position` reports a nucleotide offset (0-based), residue position
(1-based), the reference residue, the start/body/stop role, and overlapping
reference peptide **observations**, including their review quality. Its
`reference_ms_support` flag excludes low-quality and mixed reviews. It neither
validates REF nor calculates a mutant sequence; use the patient's genome and
RNA evidence for those steps. An indel or splice event can affect residues far
beyond its genomic anchor, so mapping only the anchor is not a full event
consequence analysis.

For the TPST1 example, the native 39-aa model is
`MPSRRRGGSSLIPDVGYLSEVDCPWPEHFPKIILSKISV`. Its HLA peptide `KIILSKISV` covers
residues **31–39**, with 28 PSMs in the pinned author export. The TPST1–CRCP
candidate retains only native residues 1–13 and changes the continuation.
**This native peptide provides no direct evidence for either the retained
13-aa prefix or the fusion's altered tail.** Neither reference evidence nor an
RNA-supported mutant ORF establishes tumor specificity or mutant presentation.

## Source curation and limits

The build joins the following independently pinned inputs:

1. The author's MIT-licensed copy of the
   [GENCODE Phase I ORF reference](https://github.com/VanHeeschLab/deutsch_kok_et_al_2024/blob/b65ae0c09a53e9efe083d806d617581e6139bdde/raw/ncorf_list.xlsx),
   described in [Mudge et al. 2022](https://doi.org/10.1038/s41587-022-01369-0).
   Only stopped `uORF`/`uoORF` models from its two main model sheets are included
   (minimum 16 aa). These are historical transcript models, not a replacement
   for current patient transcript annotation. The source's legacy MS study lists
   are pooled at ORF level; they do not assign each peptide to a sample or study.
   Only exclusively whole-proteome MS/MS lists receive the `shotgun` label.
2. The same pinned author's [HLA mapping and count inputs](https://github.com/VanHeeschLab/deutsch_kok_et_al_2024/tree/b65ae0c09a53e9efe083d806d617581e6139bdde/raw),
   associated with [Deutsch et al. 2026](https://doi.org/10.1038/s41586-026-10459-x).
   These are the repository's data exports, not the paper's final supplementary
   curation tiers. Counts are peptide PSMs across MS runs. They can be shared by
   several ORF mappings and must not be summed as independent samples or RNA
   abundance. `mapping_count` is the author's filtered reference-mapping count,
   not uniqueness across every human normal proteome. Missing counts stay null.
3. The separately CC-BY-4.0-licensed
   [Wacholder community assessment archive, version 1](https://doi.org/10.6084/m9.figshare.30131869.v1),
   associated with the [published benchmark](https://doi.org/10.1038/s41467-025-68002-x).
   The archived review input is retained as version 1, rather than represented
   as all final article tables. Every compatible peptide interval is retained,
   treating I/L as indistinguishable and marking that equivalence. Controls are
   excluded. Each spectrum keeps its USI and individual reviewer scores.
   All reviewers scoring 4–5 yields `reviewed_support`; all below 4 yields
   `reviewed_low_quality`; disagreement across that threshold yields
   `reviewed_mixed`. The latter two do not pass evidence filters. Source ORF
   aliases are checked; otherwise attribution is `sequence_compatible` and
   does not identify the producing ORF. Mapping uniqueness remains unknown.
   The [2026 correction](https://doi.org/10.1038/s41467-026-73431-3) concerns figure
   labels, rather than source spectrum data.

4. Ensembl release 101 GRCh38 GTF and cDNA sequences, used only after exact
   entire-ORF path/protein validation. Source URLs and SHA256 pins are in the
   manifest; per-context transcript IDs retain their versions. Ensembl makes
   its generated data [available without restriction](https://www.ensembl.org/info/about/legal/disclaimer.html).

An unreviewed `reported` observation is searchable evidence with its original
provenance, not an independently validated identification. `require_reviewed`
and `require_unique` apply to the **same observation**; a review on one spectrum
does not confer uniqueness on another observation of that sequence. Conflicting
sources remain distinct. Catalogue membership does not automatically affect
Isovar's existing variant calls, protein reconstruction, or vaccine ranking.

The complete manually curated **2026 shotgun** tables are not included. The
public PeptideAtlas non-HLA export currently exposes an older search reference
and misses published positive controls. Its absence of a peptide is not used
as negative evidence. The source manifest states this gap explicitly. Paper
supplementary files with more restrictive licenses are not silently folded into
the separately licensed repository and archive extracts.
Complete contemporary shotgun curation is tracked in
[issue #466](https://github.com/openvax/isovar/issues/466).

## Rebuilding

From an Isovar checkout with the optional builder extra and `Rscript` installed:

```sh
python -m pip install '.[uorf-build]'
python -m examples.build_uorf_catalog --source-dir /tmp/isovar-uorf-sources --prepare-sources
```

Preparation retrieves about 246 MB of pinned author inputs and a small benchmark
CSV through HTTP byte ranges, rather than its 4.8 GB archive. It temporarily
extracts the author's ~825 MiB Rdata, exports only HLA mappings, removes that
scratch file, and streams the large count table without unpacking it. It
requires sufficient RAM to load the original Rdata. Both original and derived
input checksums are checked. A repeat build without `--prepare-sources` reads
only local files and requires no R or network access. Invalid source assignments
cannot be promoted: sequence-incompatible assignments go in the manifest's
explicit exclusion list; structurally invalid models fail validation.

Protein intervals retain every exact or I/L-compatible occurrence. Records are
sorted and gzip output has a fixed timestamp. The manifest records input pins,
output checksum, scope and statistics. Attribution and redistribution notices
are packaged alongside the dataset in `NOTICE.txt` and `AUTHOR-LICENSE.txt`.
