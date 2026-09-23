# RNA paths around a nominated SV

Isovar 1.20 added `isovar sv-rna` / `reconstruct_sv_rna`: event-directed RNA
reconstruction from an indexed BAM, upstream of the [supplied-fusion
workflow](fusion.md). It is the first slice of
[#305](https://github.com/openvax/isovar/issues/305) and
[#306](https://github.com/openvax/isovar/issues/306). No predicted fusion cDNA
is required. The output is **exploratory**: it is not the validated
supplied-fusion result, and it is not passed to Vaxrank's fusion input.

```sh
isovar sv-rna --bam rna.bam --input event.json --output candidates.json
```

The input JSON has `event_id`, `reference_name`, `sample_id`, `donor` and
`acceptor` (`FusionBreakpoint`), optional extra `regions`
(`[contig, start, end]`), `references` (`FusionReference`, all candidate models)
and non-empty `event_provenance` (the DNA call and its somatic-status basis).
`sv_rna_input_from_dict` decodes it for the API. Read filters, soft clips and
read-end trimming are the usual RNA options (`--min-mapping-quality`,
`--use-soft-clipped-bases`, `--read-end-profile`, ...).

For an intergenic adjacency, supply `"references": []` when there are no
applicable transcript models. Observed junction sequences and exploratory
ATG-to-stop ORFs can still be recovered. The result and exported ORFs carry
`reference_models_unavailable`; annotated splice, start and frame assessment
is unavailable. Event linkage means compatibility with the nominated DNA
adjacency, and cannot exclude ordinary splicing without annotation. Missing
models do not establish peptide novelty or translation, and do not authorize
filling unobserved RNA sequence from the genome.

Coordinates follow `fusion.md`: 0-based, interbase; a breakpoint is the retained
partner's boundary, oriented so the donor is 5'. Each path base has one
`(contig, position, strand)` placement in path orientation, or none.

## What it does

1. **Retrieve.** Search both breakpoint neighbourhoods (`--breakpoint-window`,
   default 1000 bases each side) first. Distinct windows take turns admitting
   new records, even when they overlap, so depth at the first breakpoint cannot
   consume the entire budget before the second is sampled. Recover their
   observed SA and mate records next, following links recursively and keeping
   only records of the requesting fragments. Then spend the remaining budget
   on reference exons, extra regions and their links. A planned but unvisited
   region never suppresses a linked-record lookup. An SA tag is a retrieval
   hint, never an alignment; unplaced unmapped mates and genome-wide alternative
   placements are not assessed. Record, query and path limits are reported in
   `limitations`, never silent. A limit may still truncate breakpoint evidence.
   `acquisition.searched_regions` retains the requested regional scope;
   `region_queries` and `hop_queries` record which queries were fetched and
   exhausted (`complete`), including the priority of each regional query.
2. **Observe.** Each sequenced segment (read group, QNAME, mate) becomes one
   query path. Supplementary records are joined only when the existing
   reciprocal SA validator links them; alternative/secondary placements stay
   separate observations. Bases come from actual CIGARs in original query
   coordinates. Hard-clipped bases stay unavailable. Bases that two pieces
   place differently (junction homology) stay in the sequence but unplaced.
   A CIGAR `N` gap shorter than 21 bases is a deletion, not an intron
   (STAR's `alignIntronMin`). Missing QUAL is kept unless `--require-base-qualities`.
3. **Seed.** Paths start at unannotated junctions: CIGAR `N` or supplementary
   splits between consecutive placed bases. With soft clips enabled, unaligned
   sequence exactly at a breakpoint also seeds a path. Joins at or across the
   event (breakpoint, event-compatible, clips) may rest on one fragment.
   Splice-ambiguous and regional joins are queued and seeded by support, if the
   path budget allows. They need `--min-alternative-fragments` (default 2),
   counted over every read that makes the join. Placements of the adjacency
   (other junction-base assignments, and near-misses within
   `--max-breakpoint-shift`) are ordered by support. The strongest seeds;
   another seeds only if it is itself supported (next item). A base variant
   of one placement is pruned only against a supported best. Unannotated
   ordinary-splice or regional joins within `--annotated-junction-tolerance`
   (5) bases of an annotated junction are aligner wobble and seed nothing.
4. **Extend.** Paths grow 5' and 3' through exact overlaps (`--min-overlap`,
   default 30) that share a placed base with the path and never conflict in
   placement. There is no sequence-only joining across repeats and no bridging of
   unobserved mate inserts. At a disagreement, an alternative is pruned when
   another has at least `--min-alternative-fragments` fragments and it has
   less than `--min-alternative-fraction` (0.1) of the best. A base or
   small-indel alternative (placed within 20 bases of the best, or unplaced)
   also needs `--min-local-variant-fraction` (0.5) of the best. Systematic
   long-read errors, often 20–40% of reads, therefore do not fork paths, while
   a heterozygous allele does. Other alternatives (splice choices) fork. Pruned
   branches are listed. Once most of a step's reads have ended, the path end is
   re-queried, so one long read cannot decide the rest of a path alone.
   `--no-assembly` uses only reads that span the seed junction.
   With assembly enabled, both regional and seed-spanning reconstructions are
   considered; a seed-spanning path contained in a regional path is collapsed
   into it. This preserves witnessed flanks when competing regional extensions
   block their assembly (for example at a repeated genomic placement), without
   suppressing a regional isoform supported by shorter reads. Each path lists
   its `reconstruction_scopes` (`regional`, `seed_spanning`); containment adds
   a scope without counting its reads again. The same support thresholds and
   global path budget apply to both scopes. A retained path can still be a
   consensus: full-ORF witnesses remain a separate evidence measure.
5. **Translate.** A frame is transferred from each exact, collinear CDS match
   of at least `--min-anchor-bases` (18) and read downstream through observed
   sequence to the first stop or path end. The continuation need not match any
   annotated ORF. A substitution keeps the frame. A reading is reported only when
   it leaves its model at a reported junction or breakpoint clip. Readings
   leaving at an indel, a sequencing error or another model's splice junction
   are counted in `readings_departing_elsewhere`. Identical proteins from
   several models are grouped with every model's frame evidence.
6. **Explore other starts.** Each path also has `exploratory_orfs`, separate
   from `translations` and `frame_status`. It enumerates observed ATGs whose
   potential translation crosses an event-related junction, including
   read-through-ambiguous joins and breakpoint clips. The default minimum is
   15 aa (`--min-orf-amino-acids`); at most 100 candidates per path are retained
   in start-offset order (`--max-orf-candidates`). `candidate_limit_reached`
   reports truncation. Non-ATG initiation and starts outside the retained RNA
   are not searched. A path end without a stop remains partial.

   The ATG-to-stop interval may cross either boundary of unplaced junction
   sequence: `donor_to_unplaced` or `unplaced_to_acceptor`, or both. Adjacent
   placed flanks are `flank_to_flank`. `junction_crossings` records these
   classes and `termination_only` when only the stop codon crosses a boundary;
   such a candidate need not contain a new amino-acid sequence. Codons may
   straddle boundaries. ORFs wholly within one flank or the unplaced interval
   are excluded. Unplaced bases can be insertion, homology or clipped sequence;
   these classes do not verify a DNA insertion. Translation continues to use
   the [standard genetic code](https://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi#SG1).

   `start_context` reports the observed flanks, including the −3 and +4 bases
   when present. This is not a Kozak score or an initiation prediction. Exact
   ATG placements are compared with every supplied reference model and labeled
   annotated start, 5′ UTR, internal CDS, 3′ UTR or noncoding transcript.
   `annotated_start` means at least one supplied model annotates that codon;
   no matching model does not establish that a start is globally unannotated.
   `reference_comparisons` translates from the same reference ATG and reports
   its shared amino-acid prefix. For a partial RNA candidate, a difference
   from the complete reference ORF may only reflect truncation. Global peptide
   novelty and protein expression are not established.

   `start_evidence` adds transcript-specific origin tiers: exact annotated CDS
   ATG (priority 1), 5′-UTR ATG (2), or sequence-only start (4). The latter
   retains intronic, antisense, noncoding, junction-created and unplaced
   subtypes. Priority 3 is reserved for qualified splice-linked intronic
   inclusion and is not yet emitted. Annotation never proves initiation.
   Overlapping models with different tiers remain ambiguous, with null overall
   priority; each model's assessment is retained. For example, TPST1's 30-aa
   candidate is a 5′-UTR start in ENST00000304842 but also overlaps supplied
   noncoding isoforms, so it has no single isoform-independent priority.
   `annotate_orf_start` exposes this annotation independently of reconstruction;
   `summarize_orf_start_evidence` applies the same conservative rule across
   occurrences. See the [tier definitions and splice-inference plan](orf-start-evidence-plan.md).

   Each candidate's `full_interval_support` requires an exact complete
   nucleotide witness, including the stop when present, from a built
   observation that makes a crossed junction. Placements must be compatible,
   with shared anchors on both available sides. If an ORF starts or stops
   inside unplaced junction sequence, the same observation must also match
   the sequence and compatible placements between that ORF and the other
   available flank. `witnesses[].junction_links` records each supported
   junction's path query interval and this potentially larger observation
   interval, separately from the ORF's own interval. Alignment homology may
   require looking beyond the first flank base for a shared anchor. Mismatches
   in this extra linkage interval conservatively exclude a witness, even if
   its ORF bases match. No pair of partial reads is
   counted as one full witness. Counts distinguish segments, RG/QNAME fragments
   and RG/cell/UMI labels within this input source; these labels are not proof
   of independent molecules. Missing tags yield null molecule counts. Witness
   query intervals and observation IDs link to the original records. Counts
   cover built direct observations, not every unbuilt read beyond resource
   limits. A path stop with zero witnesses is still only assembled sequence.

## Filterable read evidence

`record_evidence` is keyed by the same IDs as `original_records` and separates
MAPQ, QUAL availability, query/aligned-query lengths and selected native tags.
`alignment_metadata` retains RG and PG header records once per result, including
platform, library and basecaller/model provenance when supplied by the producer.
Absent tags and MAPQ 255 are null; zero remains zero. Missing QUAL stays
optional (`--require-base-qualities` opts out). Original SAM retains every tag
and the actual qualities, including bases outside the candidate interval.

| Field | Meaning and use |
|---|---|
| `mapping_quality` | Confidence in alignment placement; existing `--min-mapping-quality` filter |
| `mapping_quality_raw`, `sam_flag` | Original MAPQ (including 255) and SAM flags |
| `base_qualities_available` | Whether per-base Phred qualities are stored; not the same as MAPQ |
| `original_base_qualities_tag_present` | OQ exists; it is not automatically substituted for QUAL |
| `qs`, `dx`, `pi`, `sp` | ONT mean read Q, duplex status, split-read parent and signal offset; dx is 1 for duplex, 0 for unpaired simplex, -1 for a simplex parent of duplex |
| `NM`, `mg` | Edit distance and gap-compressed alignment identity (%); true variants contribute to differences |
| `AS`, `NH`, `HI`, `nM` | Alignment score, alignment multiplicity/index and STAR mismatch count |
| `ms`, `s1`, `s2`, `dv`, `de`, `rl`, `tp` | minimap2 alignment/chaining scores, divergence, repetitive-seed length and alignment type |
| `rm` | pbmm2 trimmed overlapping query matches between alignments |
| `rq`, `np`, `ec` | Predicted read accuracy, complete insert passes, effective subread coverage when present |
| `ic`, `is`, `im` | Iso-Seq consensus input count, associated read count, input names; not independent molecule counts |
| `CB`, `UB`, `XM`, `RG` | Cell, UMI and read-group provenance; XM is raw after Iso-Seq tag and corrected after correct |
| `rc` | Iso-Seq predicted real-cell flag; neither a malignancy label nor base accuracy |
| `ff` | PacBio CCS failure bit mask, when available |
| `CR`, `CY`, `UR`, `UY`, `RX`, `QX`, `MI`, `PG` | Raw cell/UMI bases and qualities, SAM molecular barcode/quality, molecule identifier and producing program |

These are native producer-specific fields, not interchangeable quality scales.
For example, CCS may use `rq=-1` for unpolished reads. `XM` has an Iso-Seq meaning
only in an Iso-Seq context. The common evidence extractor is also available as
`isovar.read_metadata.record_evidence(read)` for ordinary pysam records. It reads
selected tags only, avoiding signal, kinetics and modification arrays.

Small variants and SVs use `ReadCollector.alignment_filter_reason`: unmapped,
vendor-QC-failed, disallowed duplicate/secondary, missing sequence/name, low
MAPQ and (when requested) missing QUAL records are excluded before assembly.
The existing numeric MAPQ filter retains 255, which STAR uses for unique
mappings. Exported `mapping_quality=null` correctly avoids claiming Q255;
consult the raw value, NH and producer metadata when choosing another policy.

Python callers can apply the same additional policy in either pipeline:

```python
from isovar import ReadCollector

# An example alignment-multiplicity policy: explicitly retain missing NH.
collector = ReadCollector(
    read_filter=lambda read: not read.has_tag("NH") or read.get_tag("NH") == 1)
```

The predicate receives each eligible original record before mate merging or SV
segment construction. Exceptions propagate. SV output records whether a custom
filter was used; record its configuration in `event_provenance` for reproduction.
Default collection does not call `record_evidence` or impose new tag thresholds.
Filters can similarly inspect `rq` or `qs` with a stated missing-value policy and
known producer. They must not infer a base-specific Phred score from these tags.

For example, downstream code can inspect each full-interval witness's
`observations[id].records`, join those IDs to `record_evidence`, and require a
chosen MAPQ or `tags.mg` threshold. Missing values need an explicit retention
policy. Filtering witnesses requires recounting distinct fragments and labels;
the original aggregate count does not survive arbitrary filtering. No combined
weight, synthetic base quality or platform-specific tag threshold is imposed.
`missing_quality_segments` counts witnesses containing any record without QUAL;
it is not a per-base quality calculation for the candidate interval.

Definitions: [SAM](https://samtools.github.io/hts-specs/SAMv1.pdf),
[Dorado](https://software-docs.nanoporetech.com/dorado/latest/basecaller/sam_spec/),
[Dorado duplex](https://software-docs.nanoporetech.com/dorado/latest/basecaller/duplex/),
[CCS filtering](https://ccs.how/faq/reads-bam.html),
[minimap2](https://github.com/lh3/minimap2/blob/master/minimap2.1),
[PacBio BAM](https://pacbiofileformats.readthedocs.io/en/13.1/BAM.html),
[Iso-Seq](https://isoseq.how/isoseq-tags.html), and
[pbmm2](https://github.com/PacificBiosciences/pbmm2).
Initiation context: [Kozak's mutagenesis study](https://pubmed.ncbi.nlm.nih.gov/3943125/).
The [Sid neo-ORF report](../figures/osteosarc/neo-orfs-and-long-reads.md)
compares candidate sequence support with what is known about initiation.

## Independent evidence axes

- **Sequence:** each junction's `direct_segments` and `direct_fragments`
  count segments/templates whose **own** alignment makes that join,
  however their other bases compare with the assembled consensus. That keeps
  noisy long reads. Every assignment of the nominated adjacency's junction
  bases supports its breakpoint junction; `direct_junction_sequences` lists
  them, including CIGAR insertions from reads beyond the segment build cap.
  `direct_cell_umi_support` separately reports input/sample/library-scoped
  cell/UMI labels and unresolved evidence; these are not independent molecules.
  `direct_molecules` is a deprecated label-count alias, `null` if any supporting
  segment lacks a resolved label or declared library scope. A read whose build failed (e.g. ambiguous bases) is
  not counted. `linked_interval` marks path bases co-observed with a novel
  junction in single reads that observe every path base in between (a read
  skipping or adding an exon contributes nothing). Intervening unplaced bases
  must agree in length and sequence too. `sequence_evidence` counts the segments whose
  overlaps voted for the path. Wider context is an assembly hypothesis, not
  proven phase. Mates, supplementary pieces and duplicate records of one
  fragment count once; a QNAME in another read group is another fragment.
- **Frame:** `frame_status` is `translated` (one protein, no competing model),
  `ambiguous` (several proteins, or a noncoding model sharing the frame
  anchor), `reference_protein_only`, `unresolved` (only noncoding or UTR
  anchors) or `no_gene_anchor`. `complete_5prime` is set only when the
  annotated start codon is observed; otherwise the upstream frame is assumed.
  `translation_observed` is always false: RNA is not protein evidence.
- **Event linkage:** every junction has a `relation`. Values, strongest first:
  - `breakpoint_junction`: the join reproduces the adjacency. Unplaced junction
    bases may belong to either partner, and an aligner may place up to
    `--max-breakpoint-shift` (10) homologous bases past one breakpoint (a
    negative `breakpoint_assignment`; `0 ≤ d + a ≤` unplaced bases). The
    homology is not checked against a genome. A join just off that diagonal,
    with one side past a breakpoint, is a noisy event join, not a regional one.
  - `event_compatible_junction`: donor side to acceptor side, but not at the
    breakpoint, as when splicing removes an intronic breakpoint.
  - `breakpoint_clip_partner_unplaced`: unaligned sequence at a breakpoint.
  - `splice_ambiguous_event_junction`: crosses the event, but ordinary forward
    splicing or read-through of the reference could make it too. That means
    any such non-breakpoint join, and every placement of an adjacency whose
    diagonal includes an annotated exon end joined to an annotated exon start
    downstream on the same strand (read-through of adjacent genes).
    `breakpoint_assignment` still shows the match.
  - `regional_novel_junction`: an unannotated join near the event which
    does not cross it.

  Annotated junctions which cross the event are listed in
  `competing_annotated_junctions` and never become candidates.
  `somatic_causation_proven` is always false.

Top-level `status` is `event_linked_candidates` (one of the first three
relations), `splice_ambiguous_candidates`, `regional_candidates_only` or
`no_candidate_paths`. The output schema is `isovar.sv_rna_candidates.v2`
(1.21): v1's `spanning_*` fields became `direct_*`, and whole-path
containment counts became voting counts. A missing path
is not evidence against the event. Peptides are windows absent from the
supplied reference proteins only; this is not proteome novelty, presentation
or immunogenicity. Results retain the SAM text of cited records, excluded-record
reasons, segment-path notes and effective parameters.

## Read lineage

`junctions[].direct_read_lineage` and exploratory ORF
`full_interval_support.read_lineage` report Dorado sequencing signal ancestry
among the respective supporting segments. `segment_ids` makes the denominator
explicit; `resolved_signal_groups` counts only resolved groups, and
`unresolved_segments`/`status_counts` expose missing or conflicting evidence.
`read_lineage.segments` at the top level records each consulted segment's
producer, status, `pi` parent, `dx` class, `sp` offset and resolved signal group.

Explicit Dorado split children can share a signal group; duplex parents and
consensuses remain unresolved because `dx` does not provide their correspondence.
Unknown producers also remain unresolved. This does not change segment/fragment
counts, reconstruction thresholds, alignment compatibility, or full-ORF witness
requirements. Signal groups are **not independent molecule counts**. They are
scoped to one input and one read group; existing barcode/UMI label fields remain
separate. See the [lineage contract and primary sources](ont-read-lineage.md).

## Cell/UMI evidence

Junction `direct_cell_umi_support` and exploratory ORF
`full_interval_support.cell_umi_support` use the same aggregation policy on
their own supporting segments. The former includes direct junction support;
the latter requires exact full-interval witnesses. A read which spans only a
junction cannot supply full-ORF evidence.

`observed_labels` counts distinct reported label pairs within their declared
scope. `segment_ids`, `unresolved_segments`, `unknown_library_segments` and
`status_counts` expose the denominator and its limitations. `complete_label_count`
is `null` unless every witness has a label and known library scope. The existing
`direct_molecules` and `full_interval_support.molecule_labels` scalars alias
this complete label count. `independent_molecules` remains `null`: label
collisions and producer-specific clustering cannot be resolved by string equality.
No reconstruction threshold, sequence vote, raw fragment count or ONT signal
group count uses these labels.

The top-level `cell_umi_evidence` contains the versioned policy and each consulted
segment's scope, label, contributing UMI tags, XM semantics and resolution status.
Visible mates are consulted for conflicts even when only one mate supports the
reported interval. See [the scope and tag contract](cell-umi-evidence.md).

## Scale

Per-base observations are built only when needed: for segments whose CIGAR
(including inserted bases) or split structure crosses the event, then for lazily seeded joins in support
order, then for reads that extend a path end. Each path end builds at most
`--max-extension-segments` (200). Reads covering most of the end's window
(e.g. both sides of a junction) come first, since an extender must overlap it
all, then those reaching furthest. Reaching the cap is reported as
`extension_segment_limit`. A queued join builds at most
that many of its reads (`seed_segment_limit`). Junction counts use a cheap
CIGAR index and need no building. Placement tuples are shared across reads.
`observation_counts` reports eligible and built segments. On Sid ONT T1
TPST1–CRCP (8k records), 1.20.0 took 40 s and 1.9 GB and hit the path cap
while counting 0 of 12 noisy FOXO3 reads. 1.21 takes 5 s and 0.3 GB. The
65k-record ONT ATP5MG–KMT2A regions (previously ~15 GB, out of reach) take
about 45 s and 1 GB.

## Regression data and limits

Synthetic tests cover both strands, spliced and exact breakpoints, observed
homology, SA-only claims, distant mates, secondary placements, read-group
collisions, missing QUAL, assembly on/off, deep noisy coverage, forks,
limits and noncoding continuations. Real records:

- **ATP5MG--KMT2A** (Sid, STAR `N` join): reproduces the validated
  `MAQFVRNLVEKTPALVNG*` translation and its noncoding ambiguity. With soft
  clips, a read-through TruSeq adapter forms a second path until the explicit
  adapter profile trims it.
- **BCR--ABL1** (K562 Iso-Seq): the fixture lacks the SA partners, so no
  placed junction is invented. With soft clips, each CCS read yields a
  breakpoint-clip path whose donor-frame peptides include the validated
  junction peptides.
- **TPST1--CRCP, FOXO3, PARD3B** (Sid ONT split reads): the breakpoint junction
  is recovered, including PARD3B's 12-nt insertion. As in the supplied
  analysis, no coding frame crosses it.
- **Sid T1 long reads** ([fixtures](../tests/data/fusions/long-read/README.md)):
  PacBio TPST1–CRCP reads with the 8-nt homology placed on CRCP are the
  breakpoint junction (20 fragments). Noisy ONT FOXO3 reads directly support
  their junction (12 fragments, 8 observed CB/UB labels with unknown library scope).
  The PacBio fixture retains CB/XM, but its disconnected program history does not
  establish corrected XM labels. Neither fixture supplies a complete molecular
  denominator. The PacBio ATP5MG–KMT2A join is
  read-through-ambiguous.

The ATP5MG–KMT2A join appears in every Sid long-read product, always between
annotated splice sites of adjacent same-strand genes. RNA alone does not
distinguish that from a rearrangement. Noisy long reads are not
error-corrected, and isoform alternatives fork paths up to `--max-paths`. Out of scope here:
whole-genome clip realignment, consensus calling, complex multi-breakend
events, automatic adapter-profile inference
([#309](https://github.com/openvax/isovar/issues/309)), translation-initiation
inference, germline/co-somatic attribution
([#297](https://github.com/openvax/isovar/issues/297)) and Vaxrank ranking.

## Export exploratory ORF sequences

Add `--orf-output-prefix sample.orfs` to `isovar sv-rna` to write
`sample.orfs.json`, `.tsv`, `.protein.fasta` and `.nucleotide.fasta` alongside
the required native `--output` reconstruction. The prefix must not collide
with the native output. The API equivalents are
`isovar.export_sv_rna_orfs(result)` and `isovar.write_sv_rna_orfs(export, prefix)`.
The adapter accepts `isovar.sv_rna_candidates.v2` and emits
`isovar.sv_rna_orfs.v2`. This is the exploratory SV ATG-ORF portion of
[#324](https://github.com/openvax/isovar/issues/324); annotated-frame translations
and ordinary SNV/indel exchange remain separate work.

Every exploratory hypothesis is retained, including partial ORFs and candidates
with no exact full-interval witness. Nucleotide FASTA includes an observed stop
codon; protein FASTA excludes the stop. Identical nucleotide, amino-acid and
stop-completeness tuples within an event/reference are grouped across paths.
Synonymous nucleotide alternatives remain separate. Path occurrences retain
start context, reference comparisons, frame status, junction boundary classes
and exact witness intervals. ORF and witness intervals are zero-based half-open;
junction query intervals instead identify the two flanking base offsets.

Each occurrence preserves `start_evidence`. `start_evidence_summary` agrees
with all supplied transcript/path tiers or reports `ambiguous`; missing
metadata in an older reconstruction reports `unavailable`. Neither case is
assigned a numeric priority. TSV includes `start_tier_status`, `start_tier`,
`start_priority`, and per-occurrence `start_evidence_json`. Lower numbers are
only an annotation prior: compare RNA support and quality separately, and
inspect ambiguous alternatives instead of sorting null priorities as evidence.
Tier metadata does not enter sequence IDs or evidence counts.

Candidate IDs hash `[event_id, reference_name, nucleotide_sequence, amino_acids,
ends_with_stop_codon]`, excluding sample/source. Sequence IDs hash only their
sequence. All hashes use full SHA-256 of JSON `[domain, value]`, with sorted keys,
ASCII escaping and no whitespace. Domains/prefixes are `sv_orf`, `nt`, `aa`,
`segment`, `fragment`, `signal` and `evidence`. Segment and fragment values are
`[[sample_id, source], identity]`, where identities are respectively
`[RG, QNAME, segment_bits]` and `[RG, QNAME]`; signal values use the resolved
signal-group identity in the same input scope. Evidence-set IDs hash sorted
segment IDs. Event identity is excluded from RNA member IDs so shared reads
across events can be detected. Different set IDs do not imply disjoint evidence:
compare member IDs. Different source aliases are not evidence of independence;
these hashes are not anonymization guarantees.

Support is recomputed from the union of original segment identities, never by
adding path counts. Cell/UMI and signal-ancestry summaries use the same ledgers
as reconstruction. Known label/signal counts describe resolved subsets;
`complete_label_count` is null unless every witness is labeled with known library
scope. `independent_molecules` is always null. TSV renders unknown numbers as
blank. Source/library metadata do not establish independent biological samples.
Raw SAM and read names are not copied into the evidence summaries; witness
observation IDs resolve against the native reconstruction.

Flags describe unannotated starts, unresolved annotated path frames, partial
ORFs, absent/single witnesses, missing qualities, unknown library/lineage,
shared signal ancestry, unplaced bases, stop-only crossings, uncertain event
linkage and reconstruction/search limits. Nondefault branch thresholds are
flagged and exact parameters retained. Fragment orientation is partitioned into
original-query, reverse-complement and mixed, relative to the processed input
sequence supplied to the aligner, **not established biological RNA strand**.
Preprocessing may have reversed or normalized this sequence already. Initiation, translation, peptide novelty
and interval-specific base quality remain unassessed by this adapter. No ranking,
abundance or presentation inference is performed. Sequence translation follows
NCBI standard code 1; an ATG-to-stop sequence is a hypothesis, not proof of use.

Sources: [SAM](https://samtools.github.io/hts-specs/SAMv1.pdf),
[SA tag](https://samtools.github.io/hts-specs/SAMtags.pdf),
[fusion reconstruction assessment](https://pmc.ncbi.nlm.nih.gov/articles/PMC6802306/)
and [NCBI translation tables](https://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi).

## Compare predictions with RNA protein hypotheses

```sh
isovar sv-rna --bam rna.bam --input event.json --output candidates.json \
  --predictions predictions.json
```

The prediction document supplies matching `event_id` and `reference_name`,
optional matching `sample_id`, nonempty `source` provenance (annotator, version,
annotation release and source calls), and a `predictions` list. Each entry has a
unique `prediction_id`, `amino_acids` (or `null` for unresolved), and optionally
`complete: true`. Additional effect-class and transcript metadata is retained.
For example:

```json
{"event_id":"event-1","reference_name":"GRCh38",
 "source":{"annotator":"varcode","version":"9.4.2","ensembl_release":95},
 "predictions":[{"prediction_id":"model-1","amino_acids":"MAK","complete":true}]}
```

`prediction_comparison` retains every event-linked annotated-frame translation
and exploratory ORF. It distinguishes complete candidates from partial readings
and exact fragment matches from whole amino-acid sequence equality. A shorter
complete ORF is not relabeled a matching fragment of a longer predicted protein.
`no_exact_sequence_match` reports a sequence comparison, not a verdict about
which isoform or sequencing hypothesis is correct. Unknown predicted proteins
stay unresolved; purely regional junctions cannot validate the nominated event.

Alternative path contexts sharing the same nucleotide/protein hypothesis keep
all path IDs. Full-interval witnesses are deduplicated by original segment and
query interval within the input source. Junction support is never substituted
for full-protein support. RNA alone does not establish initiation, protein
expression, somatic causation or antigen presentation.

Exploratory ORF comparison reuses the shared ORF exporter, including stable
candidate IDs, source-scoped witness unions, orientation and uncertainty flags.
Annotated-frame translations retain their frame evidence separately; junction
support is never reported as complete protein support.

## Reproduce the Osteosarc comparison

From a checkout with Varcode >=9.4.2 and Ensembl 95 installed:

```sh
python -m examples.osteosarc_sv_validation output-directory
```

This reconstructs TPST1–CRCP, PARD3B–CDKN2B-AS1/CDKN2B and FOXO3–STRADA/CCDC47:
the three strongest named rearrangement leads in the audited catalogue (393,
48 and 31 reported RNA split reads, respectively, before library/record
reconciliation). It uses selected original PacBio/ONT records, retains their
actual support counts, and reannotates the source T1 DNA calls with Varcode.
These catalogue screening counts are not the fixture's witness counts.
The runner examines both event orientations. TPST1 and FOXO3 use the pinned
Ensembl 87 models; PARD3B reuses the released three-fusion audit with Ensembl
115 models and original ONT T1 records. DNA predictions use Ensembl 95, with
all identities and reconstruction parameters preserved in the output. Gene-pair labels do
not assert a coding gene fusion.

For a full indexed **T1** RNA alignment, supply `--bam PATH_OR_URL`; `--source URL`
records the original identity when the path is a local acquisition. Reconstruction
queries the nominated regions and reference exons with the ordinary Isovar
limits. This does not imply a whole-BAM scan or exhaustive unconstrained assembly.
Inspect each report's `acquisition` and `limitations` before interpreting absence.

The output has one complete reconstruction/comparison JSON per event/orientation and a
`protein_comparison.csv`. Source records, alternative ORFs, frame assumptions,
full-interval witnesses and unavailable predicted proteins remain inspectable.

Fixture construction is saved with the tests:

- `tests/data/fusions/build_long_read.py` defines indexed acquisition and original
  long-read selection; `isovar.sid_data` records the acquisition/checksum recipe.
- `tests/data/fusions/build_osteosarc.py` builds the reference/junction corpus;
  `build_three_fusions.py` and `audit_three_fusions.py` supply the released
  PARD3B original-read fixture and reconstruction recipe.
- `tests/data/fusions/validation-events.json` pins DNA VCF URLs, SHA256 and record
  IDs. `python -m tests.data.fusions.build_validation_events --snapshot NAME
  --output events.json` regenerates it through the shared Osteosarc cache;
  `--source-directory DIR` instead verifies original local VCFs offline.
- Synthetic comparison inputs and the original-read #337 regression are in
  `tests/test_sv_rna_comparison.py` and `tests/test_sv_rna_orfs.py`.

The source distribution includes this runner, the small fusion data, and their
construction/audit scripts; the wheel contains the runtime comparison API and
shared packaged read bundle. The archive regression reconstructs the original
PARD3B candidate offline after building and extracting the actual sdist.

## ORF warning names and v1 migration (#353)

The v2 export uses these names in JSON and TSV:

| v1 flag | v2 flag | Meaning |
| --- | --- | --- |
| `reverse_complement_query_witnesses_only` | `reverse_complement_support_only` | All complete supporting fragments use reverse-complemented processed input-read observations. Never emitted for zero support. |
| `rna_polarity_unresolved` | `rna_strand_unresolved` | Biological transcription direction has not been established. Emitted for every exploratory candidate, including original-read-only and unsupported candidates. |
| `mixed_query_orientations` | `mixed_read_orientations` | Support includes both processed-read orientations, either across fragments or within one fragment. |

These are interpretation flags, not quality failures or automatic rejection
criteria. Reverse-complement support is not inherently defective in ONT or cDNA
libraries. Genomic mapping strand (SAM FLAG 0x10), processed-read orientation and
biological RNA strand are different quantities. This exporter does not infer
protocol-aware RNA strand. Orientation counts, source scoping and template
deduplication are unchanged.

Schema v2 is explicit in JSON and in the TSV `schema` column. New exports emit
only the new names. To read stored v1 exports before filtering warnings, use:

```python
import json
from isovar import normalize_sv_rna_orf_export, write_sv_rna_orfs

with open("stored.orfs.json") as handle:
    export = normalize_sv_rna_orf_export(json.load(handle))
write_sv_rna_orfs(export, "migrated.orfs")
```

The public normalizer accepts v1/v2, copies without mutating input, maps legacy
flags once, preserves unknown flags, and leaves sequence/evidence IDs and all
counts untouched. The writer uses the same migration and always writes v2.
Stored v1 TSV warning columns can use the mapping above; regenerate TSV from
the accompanying JSON when available. Consumers supporting both versions must
normalize before filtering, rather than silently looking only for old names.

See [pysam input sequence orientation](https://pysam.readthedocs.io/en/latest/api.html#pysam.AlignedSegment.get_forward_sequence)
and [ONT adapter-based orientation](https://epi2me.nanoporetech.com/workflows/wf-single-cell/wf-single-cell-report.html).
