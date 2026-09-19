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

Coordinates follow `fusion.md`: 0-based, interbase; a breakpoint is the retained
partner's boundary, oriented so the donor is 5'. Each path base has one
`(contig, position, strand)` placement in path orientation, or none.

## What it does

1. **Retrieve.** Search both breakpoint neighbourhoods (`--breakpoint-window`,
   default 1000 bases each side), every exon of every reference model and any
   extra regions. Then fetch observed SA and mate locations, keeping only records
   of the requesting fragments. An SA tag is a retrieval hint, never an
   alignment; unplaced unmapped mates and genome-wide alternative placements
   are not assessed. Record, query and path limits are reported in
   `limitations`, never silent.
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
5. **Translate.** A frame is transferred from each exact, collinear CDS match
   of at least `--min-anchor-bases` (18) and read downstream through observed
   sequence to the first stop or path end. The continuation need not match any
   annotated ORF. A substitution keeps the frame. A reading is reported only when
   it leaves its model at a reported junction or breakpoint clip. Readings
   leaving at an indel, a sequencing error or another model's splice junction
   are counted in `readings_departing_elsewhere`. Identical proteins from
   several models are grouped with every model's frame evidence.

## Independent evidence axes

- **Sequence:** each junction's `direct_segments`, `direct_fragments` and
  `direct_molecules` count reads whose **own** alignment makes that join,
  however their other bases compare with the assembled consensus. That keeps
  noisy long reads. Every assignment of the nominated adjacency's junction
  bases supports its breakpoint junction; `direct_junction_sequences` lists
  them. Molecules are distinct cell barcode + UMI (`CB` with `UB` or `XM`),
  `null` when untagged. A read whose build failed (e.g. ambiguous bases) is
  not counted. `linked_interval` marks path bases co-observed with a novel
  junction in single reads that observe every path base in between (a read
  skipping or adding an exon contributes nothing). `sequence_evidence` counts the segments whose
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
  breakpoint junction (20 molecules). Noisy ONT FOXO3 reads directly support
  their junction (12 reads, 8 molecules). The PacBio ATP5MG–KMT2A join is
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

Sources: [SAM](https://samtools.github.io/hts-specs/SAMv1.pdf),
[SA tag](https://samtools.github.io/hts-specs/SAMtags.pdf),
[fusion reconstruction assessment](https://pmc.ncbi.nlm.nih.gov/articles/PMC6802306/)
and [NCBI translation tables](https://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi).
