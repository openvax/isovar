# RNA paths around a nominated SV

Isovar 1.20 adds `isovar sv-rna` / `reconstruct_sv_rna`: event-directed RNA
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
   Missing QUAL is kept unless `--require-base-qualities`.
3. **Seed.** Paths start at unannotated junctions: CIGAR `N` or supplementary
   splits between consecutive placed bases. With soft clips enabled, unaligned
   sequence exactly at a breakpoint also seeds a path. Seeds tied to the event
   are always used. Splice-ambiguous and regional junctions need
   `--min-alternative-fragments` (default 2).
4. **Extend.** Paths grow 5' and 3' through exact overlaps (`--min-overlap`,
   default 30) that share a placed base with the path and never conflict in
   placement. There is no sequence-only joining across repeats and no bridging of
   unobserved mate inserts. At a disagreement, an alternative is pruned only
   when another has at least `--min-alternative-fragments` fragments and it
   has less than `--min-alternative-fraction` (0.1) of the best. Otherwise each
   alternative becomes its own path. Pruned branches are listed.
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

- **Sequence:** `sequence_evidence` counts segments and fragments whose whole
  observation agrees with the path. Each junction counts the observations that
  span it, and `linked_interval` marks bases co-observed with the junction in
  single reads. Wider context is an assembly hypothesis, not proven phase.
  Mates, supplementary pieces and duplicate records of one fragment count once;
  a QNAME in another read group is another fragment.
- **Frame:** `frame_status` is `translated` (one protein, no competing model),
  `ambiguous` (several proteins, or a noncoding model sharing the frame
  anchor), `reference_protein_only`, `unresolved` (only noncoding or UTR
  anchors) or `no_gene_anchor`. `complete_5prime` is set only when the
  annotated start codon is observed; otherwise the upstream frame is assumed.
  `translation_observed` is always false: RNA is not protein evidence.
- **Event linkage:** every junction has a `relation`. Values, strongest first:
  - `breakpoint_junction`: the join reproduces the adjacency, allowing unplaced
    junction bases to be assigned to either partner (see
    `breakpoint_assignment`; this assignment is not checked against the genome).
  - `event_compatible_junction`: donor side to acceptor side, but not at the
    breakpoint, as when splicing removes an intronic breakpoint.
  - `breakpoint_clip_partner_unplaced`: unaligned sequence at a breakpoint.
  - `splice_ambiguous_event_junction`: crosses the event, but ordinary forward
    splicing or read-through of the reference could make it too.
  - `regional_novel_junction`: an unannotated join near the event which
    does not cross it.

  Annotated junctions which cross the event are listed in
  `competing_annotated_junctions` and never become candidates.
  `somatic_causation_proven` is always false.

Top-level `status` is `event_linked_candidates` (one of the first three
relations), `regional_candidates_only` or `no_candidate_paths`. A missing path
is not evidence against the event. Peptides are windows absent from the
supplied reference proteins only; this is not proteome novelty, presentation
or immunogenicity. Results retain the SAM text of junction-spanning and
whole-path records, excluded-record reasons and effective parameters.

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

Noisy long reads are not error-corrected. Exact overlaps fragment them into
many single-read branches, bounded by `--max-paths`. Out of scope here:
whole-genome clip realignment, consensus calling, complex multi-breakend
events, automatic adapter-profile inference
([#309](https://github.com/openvax/isovar/issues/309)), translation-initiation
inference, germline/co-somatic attribution
([#297](https://github.com/openvax/isovar/issues/297)) and Vaxrank ranking.

Sources: [SAM](https://samtools.github.io/hts-specs/SAMv1.pdf),
[SA tag](https://samtools.github.io/hts-specs/SAMtags.pdf),
[fusion reconstruction assessment](https://pmc.ncbi.nlm.nih.gov/articles/PMC6802306/)
and [NCBI translation tables](https://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi).
