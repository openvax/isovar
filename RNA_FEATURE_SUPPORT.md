# Common RNA feature support implementation

Implements the producer part of openvax/vaxrank#757. Isovar owns read/fragment
identity, local RNA-path observation predicates and sampling geometry. The
public contract must work for allele intervals, junctions and coding intervals
without dispatching on antigen category. Preserve existing producer exports and
legacy counts; attach explicit new measurement records and retain unavailable
values when old inputs do not establish the required feature or geometry.

Implement in order:

1. Validated feature, sampling-profile and support records; exact SE/PE
   opportunity calculation and read/fragment/depth/length-normalized values.
2. Count from existing Isovar read/path evidence with explicit complete versus
   partial predicates, and expose a public API for all source adapters.
3. Carry measurement identity, observed unit identities, model assumptions and
   unavailable reasons through source exports. No gene/transcript TPM conversion.
4. Exercise event/category equivalence, SE equality, PE joint/overlap behavior,
   complex-path uncertainty and normalization under simulation. Lint, full
   tests, CI, release; coordinate the Vaxrank dependency floor.

The mathematical and category-specific acceptance contract is recorded at
https://github.com/openvax/vaxrank/issues/757. It extends the effective fragment
start-position principle in https://doi.org/10.1186/1471-2105-12-323; geometric
normalization is an explicitly named experimental measurement model.

## Common measurement record

Every observation identifies the patient, sample, library, timepoint, source/evidence set, antigen occurrence, nucleotide path and genetic code; the exact feature and its coordinates; counting/ambiguity policy and producer version; and coverage, phase, uncertainty and availability. All antigen types use the same schema and normalization definitions. Unknown is not zero. RNA-only evidence is not proof of translation or tumor specificity.

The biological feature is an **allele-bearing base/interval, diagnostic junction/path, or defined peptide-coding interval**. A gene/transcript/protein is a carrier and provenance identity, not an automatic counting region or length denominator. Distinguish event support, selected-window support, assembly support and inferred abundance; do not silently substitute one for another. Keep counts for different features as separate records even when they belong to one antigen.

Required measurement families for both units:

| Measurement | Reads | Fragments |
| --- | --- | --- |
| Raw observed support | `C_read` | `C_fragment` |
| Whole-library denominator | `N_read` | `N_fragment` |
| Effective support length | `L_eff_read` | `L_eff_fragment` |
| Per million | `10^6 * C_read / N_read` (RPM) | `10^6 * C_fragment / N_fragment` (FPM) |
| Per effective kb | `10^3 * C_read / L_eff_read` | `10^3 * C_fragment / L_eff_fragment` |
| Per effective kb per million | `10^9 * C_read / (N_read * L_eff_read)` | `10^9 * C_fragment / (N_fragment * L_eff_fragment)` |

These last quantities are feature-support densities, not standard gene/transcript RPKM/FPKM or TPM. Gene/transcript abundance and producer expected counts remain separate typed measurements. Do not manufacture observed counts from TPM or convert read counts to fragment counts by dividing by two.

A read is a sequenced segment, not an alignment record. Supplementary/split alignments of one segment do not create extra reads. A fragment is the sequenced template: paired mates count once. Record duplicate removal and UMI policies; fragment counts do not automatically mean independent original molecules. Retain UMI/cell evidence separately when assessed.

For SE data, the same eligible units and acceptance rule must give `C_read = C_fragment`, `N_read = N_fragment`, `L_eff_read = L_eff_fragment`, and equality of corresponding normalized measurements. For PE data, a read must individually satisfy its support rule; a fragment can use either mate or joint pair evidence according to the declared rule. Joint-pair support must not invent individually supporting reads. Orphans, overlapping mates, orientation and mate conflicts have explicit policies.

Denominators describe the original whole library under a recorded eligibility policy, with matching units and read/fragment identity rules. The same denominator includes the host/viral analysis universe; no category-specific million. A regional/targeted extract does not establish a whole-library denominator. Supplied upstream library receipts can provide one, but candidate admission may not redefine it. Missing denominators leave RPM/FPM unavailable without changing raw support. A missing effective length leaves length-normalized values unavailable without changing RPM/FPM.

## Definition of effective support length

For a feature `F`, a particular RNA path/occurrence `h`, a library `s`, and counting unit `u` (read or fragment), define an explicit observation-acceptance rule `A_u`. It specifies the bases/junctions/phase required, anchor lengths, permitted alignments, quality and mapping filters, ambiguity treatment, and whether evidence is complete or partial.

The proposed normalization is:

`L_eff_u(F, h, s) = sum_g p_s,u(g) * sum_x a_u(F, h, x, g)`

- `g` is a sampling configuration: read lengths, fragment/insert length, mate/end selection, orientation, and any explicitly modeled quality/detection variables.
- `x` ranges over physically valid template start positions on the RNA path, including the available flanks, not merely positions inside the feature.
- `a_u` is the indicator or modeled probability that this configuration would satisfy **the same observation rule used to count support**. Start-position bias, if modeled, uses a declared dimensionless weighting with a fixed normalization; its scale cannot be arbitrary.
- `p_s,u` is a normalized library sampling distribution independent of the observed number of alternate-supporting units. For the read-unit model of PE data, sample a sequenced mate with its library frequency, rather than summing both mates while also normalizing by `N_read`. For fragments, evaluate the pair once.

Thus `L_eff` is an expected number of informative start positions (reported as effective bp), not the length of the deleted DNA, the complete transcript or the protein. Computing it requires local RNA geometry sufficient to enumerate the relevant starts; it does not always require a reconstructed full transcript. Unknown RNA ends must not be treated as cropped true ends.

Read and fragment effective lengths are separate. The PE fragment length distribution and actual sequenced bases matter: an unsequenced gap is not observed sequence. A declared pair-linkage predicate can support phase across that gap, but cannot claim its bases were sequenced.

For an idealized internal contiguous feature of `m` bases requiring complete observation plus `a_left` and `a_right` flanking bases, a fixed SE read length `r` gives:

`L_eff = max(0, r - (m + a_left + a_right) + 1)`

For a boundary, `m = 0` and both anchors must be positive. Transcript ends clip the start set; variable read lengths are averaged over the library distribution. A deletion junction uses the joined alternate-RNA coordinates. A geometric-only implementation must be labeled as such: it does not claim to correct unmodeled mapping, sequence or capture bias. A detection-corrected model must carry its model identity and validation. Do not use a geometry estimate as a fully corrected abundance without recording that assumption.

`L_eff = 0` means the stated rule cannot observe the feature with that library/model; normalized support is unavailable, not evidence of absent expression. Positive measured support with zero modeled opportunity is an inconsistency to resolve. Do not add an epsilon to create a score. Small opportunity estimates retain count uncertainty and cannot replace minimum-evidence requirements.

## Feature definitions by event/category

All rows below use the same `L_eff` definition above. Their explicit acceptance predicates determine the start positions; no category receives an arbitrary scaling constant.

| Event/category | Feature and accepted support | Effective-length rule / required boundary |
| --- | --- | --- |
| SNV / RNA base substitution | Alternate base plus declared anchors on the observed RNA path. Reference/other bases are separately counted. | Single-base interval (`m = 1`) with the declared anchors. |
| MNV / phased compound substitutions | The declared set of alternate bases on one supported haplotype. Complete-span and joint-pair evidence are distinguished from individual-site counts. | Enumerate starts satisfying the complete haplotype predicate; contiguous full-span mode uses its bounding interval and anchors. Never sum independent SNV counts to claim phase. |
| Short insertion / deletion / delins | Exact alternate sequence/path, including required boundaries and anchors; repeated-sequence equivalent placements are resolved or retained as uncertainty. | Use the alternate RNA feature span. A deletion is an adjacency with anchors, not the length of the removed reference segment. A delins uses the replacement interval and boundaries. |
| Long insertion | Separate insertion-boundary, insertion-internal selected coding interval, and complete-insertion observations. Internal sequence must have the required source/path attribution. | Enumerate each declared feature independently. Complete-insertion opportunity may be zero for a short-read library; boundary counts remain boundary counts. No union of disconnected observations is renamed full-insertion support. |
| Fusion / translocation / exon deletion or skipping | The novel RNA adjacency with the required left/right anchors and orientation. Split segments and explicitly accepted spanning pairs are separate evidence modes. | Enumerate starts on the joined RNA sequence using the junction rule. Partner-gene lengths and genomic separation are irrelevant. |
| Inversion / tandem duplication / complex or multibreakpoint SV | Each expressed diagnostic adjacency and the selected altered coding interval, retaining the RNA path, orientation, phase and alternate placements. | Per-feature start enumeration. A complete-path count requires a unit satisfying the complete-path predicate. An explicitly requested any-junction union uses a union of unit identities and start sets; it is labeled partial/event support. Neither summed junction counts nor their minimum becomes a complete-path count. |
| Aberrant splice / intron retention / alternative exon | The novel splice adjacency, retained-intron/exon coding interval, or boundary set relevant to the selected peptide, including required linkage. | Junction or interval opportunity as appropriate, on the expressed spliced/unspliced RNA path. Intronic genomic length is not automatically the denominator. |
| Frameshift, start/stop gain/loss, readthrough, other altered ORF consequence | Measure the causal RNA feature and the downstream selected coding interval separately. Sequence, frame, genetic code and linkage determine which interval can support the antigen. | Apply the underlying base/indel/junction rule to the cause and the interval rule to the selected coding sequence. Do not use novel-tail length to normalize causal-site counts or infer translation from RNA. A termination-only change does not create an unobserved peptide. |
| CTA | A specified peptide/vaccine-window coding interval on compatible transcript paths. Isoform-discriminating junctions are additional explicit features when needed. | Enumerate interval-spanning starts (plus declared attribution anchors), not all gene/transcript reads. Use the same declared window/coverage policy as for other non-mutant proteins. Gene-TPM-only input retains TPM; it cannot supply observed local counts or `L_eff`. |
| Non-mutant isoform / self-protein / other admitted reference antigen | The selected coding interval and any junction needed to distinguish the relevant isoform. | Same interval/junction rules as CTA. The category label cannot select different count units or normalize by whole-protein length. |
| ERV / repetitive endogenous antigen | The selected coding interval or diagnostic junction, with locus/sequence ambiguity retained. | Same interval/junction rule plus the declared mapping/attribution model. Family-wide mapped reads do not automatically support a specific locus or coding window. |
| Viral protein | The selected coding interval on the actual viral RNA transcript/strain; relevant splice or virus-host junctions remain explicit features. | Same interval/junction rule using the viral RNA path and the same library denominator as host targets. Whole-virus reads, genomic copy number and unobserved coding regions do not substitute. Polycistronic/overlapping ORFs retain distinct coding mappings. |
| ncORF / alternative frame | The selected nucleotide coding interval plus an independently supported frame/cause where required. | Same local interval/junction opportunity; overlapping frames may share RNA support without RNA proving which frame was translated. Mutation-derived ncORFs retain their causal-event feature too. |

A peptide shared by several indistinguishable isoforms can have sequence-level RNA support while its isoform origin remains unknown. Count an observation once within that target-equivalence group. Do not add the same observation once per transcript/protein attribution. Where compatible paths give different opportunities, retain path-specific values or use an explicit mixture model; do not silently average their lengths. Merely compatible reads, directly diagnostic observations and fractional/expected allocations remain distinguishable measurements. Existing Isovar allocation/evidence identities should be reused.


## API and CLI

`RNAFeature(identifier, path_id, intervals, kind, mode="all", path_length=None)`
uses **RNA-path nucleotide coordinates**, including any required anchors.
`RNASamplingProfile` carries independently supplied read/template geometry and
optional original-library denominators. `effective_rna_feature_lengths` computes
both lengths. `normalize_rna_feature_support(record, profile)` can add or replace
independently supplied library statistics on a validated observation while
retaining its original feature, unit identities and measurement lineage.
`measure_rna_feature` accepts existing `OriginalPathMatch` objects;
`measure_rna_path_features` obtains those matches from the original BAM using
Isovar's common read policy and optional reusable RNA library index.

```sh
isovar rna-features --bam original.bam --features features.json \
  --sample-id biopsy-T2 --source-id library-A \
  --sampling-profile library.json --output support.json \
  --library-cache rna-index
```

`features.json` contains `{"paths": [...]}`. Each path has exactly `sequence`
(the supplied nucleotide path), `positions` (one `[contig, zero_based_base,
strand]` per base; inserted bases are null), and `features` (a list of the
`RNAFeature` fields above). All features on a path use one `path_id`. Source
annotations/reconstruction supply these coordinates; the feature command does
not infer a coding frame or a gene's abundance from a category label.

Example `library.json` for independently verified 2×100 bp sequencing and a
single 250 nt RNA-template geometry (a real distribution supplies multiple
weighted entries):

```json
{
  "provenance": "original-library-QC-sha256-and-count-policy",
  "geometries": [{"segments": [[0,100],[150,250]], "template_length": 250, "weight": 1}],
  "library_reads": 2000000,
  "library_fragments": 1000000,
  "denominator_scope": "whole_library"
}
```

The CLI writes `isovar.rna_feature_measurements.v1`, with a `measurements` list.
Every item is a self-contained `isovar.rna_feature_support.v1` record, validated
by `validate_rna_feature_support`: counts, identities, feature, geometry,
normalizations, unavailable reasons and provenance. A path model supplies no
observations by itself. Exact coordinate-compatible matches do not assert
unique transcript/locus attribution; known original-placement uncertainty is
retained. No source transcript/isoform ambiguity is resolved by this interface.

The uniform geometry model excludes sequence/mapping/capture biases and does
not turn local evidence into calibrated transcript abundance. When full path
ends are unknown, `path_length=None` explicitly selects an unbounded local-path
approximation. Use a complete known path length for boundary correction; never
use the edges of a cropped reconstruction as biological transcript ends.

Local matching checks only the declared required bases and their coordinate
placements. A mismatch outside a feature cannot erase an exact local witness.
Disjoint exact runs from the same original placement may support a compound
feature; a mismatch inside a required interval prevents complete support.
Alternative placements or offsets on a repeated path are never combined to
invent phase. Required reference/alternate feature bases must be resolved:
unknown model bases are an input error, not measured zero support.
