# RNA paths around a nominated SV

`isovar sv-rna` asks: *given a structural variant called from DNA, what RNA
sequences do the reads actually show around it, and could any of them encode a
new protein?* It searches an indexed RNA BAM around the two breakpoints,
assembles the RNA paths the reads support, and reports each path's junctions,
reading frames and candidate proteins, together with the evidence behind them.

Use it when you have an SV call but not a fusion transcript. If you already have
a sequence-resolved fusion transcript from another tool, validate it with
[`isovar fusion`](fusion.md) instead. Both commands produce the same
[RNA path format](#the-rna-path-format), so downstream code can read either.

The output is **exploratory**. It shows what the RNA supports; it does not prove
the SV causes the RNA, that a protein is made, or that any peptide is presented.

## Quick start

```sh
isovar sv-rna --bam rna.bam --input event.json --output candidates.json
```

`event.json` names the event and the transcript models to compare against:

```json
{
  "event_id": "GENE1--GENE2",
  "reference_name": "GRCh38",
  "sample_id": "tumor-1",
  "donor":    {"contig": "1", "position": 2500, "strand": "+"},
  "acceptor": {"contig": "2", "position": 5500, "strand": "+"},
  "regions": [["1", 2400, 2600]],
  "references": [
    {"transcript_id": "ENST00000000001", "reference_name": "GRCh38", "annotation": "Ensembl 115",
     "contig": "1", "strand": "+", "exons": [[1000, 1100], [2000, 2150]],
     "sequence": "GAGGCACGACCA...", "cds_start": 20, "cds_end": 200}
  ],
  "event_provenance": {"caller": "manta", "version": "1.6", "somatic_basis": "absent in matched normal"}
}
```

- **`donor` and `acceptor`** are the breakpoints: the last retained base boundary
  of each partner, oriented so the donor is 5′ of the join. Coordinates are
  0-based.
- **`references`** lists every transcript model that might be involved, with
  exons, spliced sequence and CDS. Include noncoding isoforms too, so they can
  compete. Use `"references": []` for an intergenic event. Junctions and
  start-codon ORFs are still found, but annotated frames and splice assessment
  are unavailable, and the result says `reference_models_unavailable`.
- **`regions`** (optional) adds search intervals. Both breakpoint neighbourhoods
  and the reference exons are always searched.
- **`event_provenance`** is required and is copied into the output verbatim:
  record where the DNA call came from and why it is considered somatic.

Unknown keys are rejected, so a misspelled field fails immediately. From Python,
`isovar.sv_rna.sv_rna_input_from_dict(data)` decodes the same document and
`isovar.reconstruct_sv_rna(bam, source=..., **inputs)` runs the reconstruction.
All the usual RNA read options apply, such as `--min-mapping-quality`,
`--use-soft-clipped-bases` and `--read-end-profile`.

## Reading the result

Start at the top-level **`status`**:

| `status` | Meaning |
|---|---|
| `event_linked_candidates` | At least one path has a junction that reproduces or is compatible with the nominated adjacency |
| `splice_ambiguous_candidates` | The best paths cross the event, but ordinary splicing or read-through could also explain them |
| `regional_candidates_only` | Only unannotated junctions near the event that do not cross it |
| `no_candidate_paths` | No path was assembled. This is not evidence against the event: check `limitations` and `acquisition` |

Then look at each entry in **`paths`**. A path is one assembled RNA sequence:

- **`sequence`** is the RNA sequence in 5′→3′ path orientation.
- **`junctions`** lists every non-reference join in it, with its `relation` to the
  event (below) and its direct read support.
- **`frame_status` and `translations`** give proteins from annotated reading
  frames carried across the junction.
- **`exploratory_orfs`** lists start-codon-to-stop ORFs that cross an
  event-related junction, in any frame.
- **`end_reasons`** says whether each end of the path is a real end of the
  observed RNA or was cut short by a limit.

### Junction relations

Each junction's `relation` says how it relates to the nominated event, strongest first:

| `relation` | Meaning |
|---|---|
| `breakpoint_junction` | The RNA join reproduces the DNA adjacency, allowing up to `--max-breakpoint-shift` (10) homologous bases to sit on either side |
| `event_compatible_junction` | Donor side joined to acceptor side, but not exactly at the breakpoint, as when splicing removes an intronic breakpoint |
| `breakpoint_clip_partner_unplaced` | Unaligned (soft-clipped) sequence starting exactly at a breakpoint; needs `--use-soft-clipped-bases` |
| `splice_ambiguous_event_junction` | Crosses the event, but splicing or read-through of the unrearranged genome could make the same join (for example adjacent genes on one strand) |
| `regional_novel_junction` | An unannotated join near the event that does not cross it |

Annotated junctions that happen to cross the event are listed in
`competing_annotated_junctions` and never become candidates.
`somatic_causation_proven` is always false: RNA cannot establish that the SV caused the join.

### How to judge a candidate

Isovar reports sequence, frame and event linkage as separate kinds of evidence.
Check each on its own:

1. **Is the junction real?** `direct_fragments` counts read fragments whose *own*
   alignment makes that join. One or two fragments is thin evidence.
   `direct_junction_sequences` shows every junction-base variant seen.
2. **Is it linked to the event?** Prefer `breakpoint_junction` and
   `event_compatible_junction`; treat `splice_ambiguous_event_junction` with care.
3. **Is the wider sequence observed, or assembled?** Bases far from the junction
   come from overlapping reads (`sequence_evidence`), not necessarily from one
   molecule. `linked_interval` marks the bases seen in the same reads as the junction.
4. **Is the protein justified?** `frame_status` below. For an exploratory ORF,
   `full_interval_support` counts reads that contain the *entire* ORF and a crossed
   junction; zero means the ORF is only assembled sequence.
5. **Was anything cut short?** Any reason in `end_reasons` other than
   `observations_end` means that end was truncated. Also check the top-level
   `limitations` (for example `record_limit`, `path_limit`).

Counts of segments, fragments, cell/UMI labels and ONT signal groups are
different units. None of them is a count of independent RNA molecules.

## The RNA path format

`isovar sv-rna` (schema `isovar.sv_rna_candidates.v4`) and `isovar fusion`
(schema `isovar.fusion_rna.v2`) share these fields. Every interval is 0-based
and half-open (`[start, end)`) in path coordinates.

| Level | Fields |
|---|---|
| Result | `schema`, `event_id`, `reference_name`, `sample_id`, `status`, `paths`, `parameters` |
| Path | `path_id`, `sequence`, `junctions`, `frame_status`, `translations` |
| Junction | `query_interval`: the unplaced bases between the two partners (empty for a direct join); `left`/`right`: the placed bases on either side as `[contig, position, strand]`; `unplaced_bases`; `relation`; `direct_fragments` |
| Translation | `translation_start`, `translation_end`, `amino_acids`, `ends_with_stop_codon`, `complete_5prime`, `transcript_ids`, `frame_evidence`, `candidate_peptides` (`sequence`, `protein_interval`), `translation_observed` (always false) |

`path["sequence"][start:end]` for a junction's `query_interval` equals its
`unplaced_bases`. Each tool adds fields of its own; see below and the
[fusion guide](fusion.md).

## Frames and translations

A frame is transferred from each exact, collinear match to an annotated CDS of
at least `--min-anchor-bases` (18) and read through the observed sequence to the
first stop or the end of the path. The protein after the junction need not match
any annotated ORF. A substitution keeps the frame. A reading is reported only
when it leaves its transcript model at a reported junction or breakpoint clip.

| `frame_status` | Meaning |
|---|---|
| `translated` | One protein, and no competing noncoding model shares the frame anchor |
| `ambiguous` | Several proteins, or a noncoding model shares the anchor: do not pick one silently |
| `reference_protein_only` | The frame reads straight through the reference; no junction changes it |
| `unresolved` | Only noncoding or UTR anchors |
| `departs_elsewhere_only` | A coding anchor exists, but every reading leaves its model at an indel, a sequencing error or another splice before a reported junction (counted in `readings_departing_elsewhere`) |
| `no_gene_anchor` | No exact coding anchor |

`complete_5prime` is true only when the annotated start codon itself is in the
RNA; otherwise the protein assumes the annotated upstream frame. Identical
proteins from several models are grouped, with every model's `frame_evidence`.
`candidate_peptides` lists windows of `--peptide-lengths` (8–11) absent from the
supplied reference proteins. That is not proteome-wide novelty.

## Exploratory ORFs

`exploratory_orfs` looks beyond annotated frames. It lists every observed ATG
whose ORF crosses an event-related junction, in start order. The defaults keep
ORFs of at least `--min-orf-amino-acids` (15) and at most
`--max-orf-candidates` (100) per path; `candidate_limit_reached` reports
truncation. Non-ATG starts are not searched. An ORF without a stop is partial.

For each candidate:

- **`junction_crossings`** says which boundary of the junction the ORF crosses:
  `donor_to_unplaced`, `unplaced_to_acceptor`, or `flank_to_flank` for a direct
  join. `termination_only` is set when only the stop codon crosses; such an ORF
  may add no new amino acids. Unplaced bases can be an insertion, homology or
  clipped sequence.
- **`start_context`** shows the bases around the ATG, including −3 and +4. This is
  not a Kozak score.
- **`reference_comparisons`** translates from the same ATG in each supplied
  model and reports the shared amino-acid prefix. A partial RNA ORF can differ
  from a complete reference ORF simply because it is truncated.
- **`start_evidence`** ranks where the start lies, per transcript: 1 for the
  annotated CDS start, 2 for a 5′-UTR ATG, 3 for an intronic ATG with a
  splice-linked inclusion witness, and 4 for anything else. Models that disagree
  leave the overall priority null. Priorities are an annotation prior, not a
  probability of translation. See [ORF start evidence](orf-start-evidence.md).
- **`full_interval_support`** counts only observations that contain the whole
  ORF, including its stop when present, and also make a crossed junction, with
  compatible placements on both sides. Two partial reads never add up to a full
  witness. `witnesses[]` lists the observation IDs and intervals, including a
  wider `junction_links` interval when the ORF starts or stops inside unplaced
  junction bases.

`isovar.annotate_orf_start` and `isovar.summarize_orf_start_evidence` expose the
start annotation outside reconstruction. The [Sid neo-ORF report](../figures/osteosarc/neo-orfs-and-long-reads.md)
compares candidate sequence support with what is known about initiation.

## Link an intronic start to an observed transcript path

An ATG inside an intron could be real if the RNA shows that intron segment
spliced into a mature transcript. Intronic `start_evidence` assessments carry a
`splice_inclusion` record for this. Starts in other regions have none and keep
`splice_inference="not_assessed"`. `isovar.annotate_orf_inclusion` runs the same
assessment on retained observations.

An intronic start reaches priority 3 only when one full-ORF witness links it to
both the event and a qualifying CIGAR `N` splice:

- **Cryptic donor or acceptor use**: the start's own aligned segment joins an
  annotated exon boundary. Two such joins on one observation suggest exonization.
- **Retained intron**: the segment spans both exon/intron boundaries, the same
  observation has another annotated splice, and other reads splice this intron out.
- **Rearranged intronic segment**: joined to an annotated partner exon that has its
  own linked annotated splice. A rearrangement alone does not prove maturation.

The linked interval must match the original read sequence and placements, with:
- at least `--inclusion-min-splice-anchor-bases` (8) bases on each side of the splice;
- MAPQ of at least `--inclusion-min-mapping-quality` (20);
- every base at `--inclusion-min-base-quality` (Q20) or higher.

The same MAPQ gate selects the competing splices used for intron retention.
Missing base qualities fail the gate, reported as unavailable rather than low. MAPQ 255 counts as a unique alignment, as
STAR writes it. For aligners that use 255 as the SAM specification's
"unavailable", pass `--inclusion-mapq-255-unavailable`; 255 then fails the gate,
and `limitations` lists `mapq_255_excluded_from_inclusion`. The gates used are
recorded in `parameters.inclusion` and in each assessment's `thresholds`.
Matched-normal comparison and splice-predictor scores are not performed and are
reported as `not_assessed`.

## How the reconstruction works

1. **Retrieve.** Search both breakpoint neighbourhoods (`--breakpoint-window`,
   1000 bases each side) first. The two windows take turns admitting records, so
   depth at one breakpoint cannot use up the budget before the other is sampled.
   Next, follow the reads' supplementary (`SA`) and mate locations, then the
   reference exons and extra regions. SA tags are only retrieval hints. Unmapped
   mates and genome-wide alternative placements are not examined. `acquisition`
   records every query and whether it finished.
2. **Build observations.** Each sequenced segment (read group, read name, mate)
   becomes one observation. Supplementary pieces are joined only when their SA
   tags reciprocally confirm the same path; secondary placements stay separate.
   Hard-clipped bases are unavailable. Bases two pieces place differently
   (junction homology) are kept but unplaced. A CIGAR `N` shorter than 21 bases is
   treated as a deletion (STAR's `alignIntronMin`).
3. **Seed.** Paths start at unannotated junctions: CIGAR `N` gaps or
   supplementary splits. With `--use-soft-clipped-bases`, clipped sequence exactly
   at a breakpoint also seeds. Joins at or across the event may rest on one
   fragment. Other joins need `--min-alternative-fragments` (2). Joins within
   `--annotated-junction-tolerance` (5) bases of an annotated junction are
   treated as aligner wobble.
4. **Extend.** Paths grow in both directions through exact overlaps of at least
   `--min-overlap` (30) bases that share a genomic placement with the path. Where
   reads disagree:
   - An alternative is pruned when the best has at least
     `--min-alternative-fragments` fragments and the alternative has under
     `--min-alternative-fraction` (0.1) of its support.
   - A base or small-indel alternative also needs `--min-local-variant-fraction`
     (0.5) of the best. This way systematic long-read errors (often 20–40% of reads)
     do not fork paths, but a heterozygous allele does.
   - Real alternatives, such as splice choices, fork into separate paths.

   Pruned branches are listed in `pruned_branches`. `--no-assembly` uses only
   reads that span the seed junction. With assembly, paths built from regional
   overlaps and from junction-spanning reads alone are both kept
   (`reconstruction_scopes`), so a witnessed flank is not lost when regional
   extension goes elsewhere.
5. **Translate and explore**, as described above.

Per-base observations are built lazily, and `--max-extension-segments` (200)
bounds the reads built to extend one path end. On Sid ONT T1, TPST1–CRCP (8,000
records) takes about 5 s and 0.3 GB, and the 65,000-record ATP5MG–KMT2A region
takes about 45 s and 1 GB.

## Options

| Option | Default | What it controls |
|---|---|---|
| `--breakpoint-window` | 1000 | Bases searched on each side of each breakpoint |
| `--max-breakpoint-shift` | 10 | Homologous junction bases allowed past a breakpoint |
| `--min-overlap` | 30 | Exact overlap needed to extend a path |
| `--min-anchor-bases` | 18 | Exact CDS match needed to transfer a frame |
| `--min-alternative-fragments` | 2 | Support needed to seed an unattributed join or prune an alternative |
| `--min-alternative-fraction` | 0.1 | Alternatives below this fraction of the best are pruned |
| `--min-local-variant-fraction` | 0.5 | Base/small-indel alternatives also need this fraction of the best |
| `--annotated-junction-tolerance` | 5 | Unannotated joins this close to annotated ones are wobble |
| `--no-assembly` | off | Only use reads spanning the seed junction |
| `--max-records`, `--max-queries`, `--max-paths` | 10000, 1000, 100 | Resource limits; reaching one is reported in `limitations` |
| `--max-extension-segments` | 200 | Reads built to extend one path end |
| `--peptide-lengths` | 8 9 10 11 | Candidate peptide lengths |
| `--min-orf-amino-acids`, `--max-orf-candidates` | 15, 100 | Exploratory ORF length and per-path cap |
| `--inclusion-min-splice-anchor-bases`, `--inclusion-min-base-quality`, `--inclusion-min-mapping-quality` | 8, 20, 20 | Splice-linked inclusion gates |
| `--inclusion-mapq-255-unavailable` | off | Treat MAPQ 255 as unavailable in the inclusion gate |

`parameters` records every threshold used, `parameters.read_collection` every
read-collection setting (including trimming and read-end profile SHA-256s), and
`parameters.inclusion` the inclusion gates. Records without a read name are
skipped and reported as `unnamed_records_skipped`.

## Export exploratory ORFs

Add `--orf-output-prefix sample.orfs` to also write `sample.orfs.json`,
`sample.orfs.tsv`, `sample.orfs.protein.fasta` and `sample.orfs.nucleotide.fasta`.
Existing files are replaced. From Python, use `isovar.export_sv_rna_orfs(result)`
and `isovar.write_sv_rna_orfs(export, prefix)`. The export (schema
`isovar.sv_rna_orfs.v3`) keeps every hypothesis, including partial ORFs and
ORFs without a full witness:

- Identical sequences found on several paths become one candidate, with each
  path listed as an occurrence. Synonymous alternatives stay separate.
- Support is recomputed from the union of original read segments, never by
  adding path counts.
- The nucleotide FASTA includes an observed stop codon; the protein FASTA does not.
- `start_evidence_summary` agrees with every occurrence's start tier, or reports
  `ambiguous` with no priority.

`uncertainty_flags` lists what to check before using a candidate:

| Flag | Meaning |
|---|---|
| `partial_orf` | No stop codon before the path ends |
| `no_full_fragment_witness`, `single_full_fragment_witness` | Zero or one read contains the whole ORF |
| `no_annotated_start_in_supplied_models` | No supplied model starts at this ATG |
| `annotated_path_frame_unresolved_or_ambiguous` | The path's annotated frame is not a single translation |
| `start_tier_ambiguous` | Transcript models disagree about the start's origin |
| `path_end_truncated` | A path end was cut short by a limit |
| `missing_base_qualities` | A witness lacks QUAL |
| `library_scope_unresolved`, `cell_umi_labels_unresolved`, `signal_lineage_unresolved`, `shared_signal_ancestry` | Label or lineage evidence is incomplete, or witnesses share an ONT signal |
| `unplaced_junction_sequence`, `termination_only_junction_crossing`, `splice_ambiguous_event_junction`, `breakpoint_clip_partner_unplaced` | The crossed junction has these properties |
| `orf_candidate_limit_reached` | The path had more ORFs than `--max-orf-candidates` |
| `rna_strand_unresolved`, `initiation_unobserved`, `translation_unobserved`, `peptide_novelty_unassessed`, `interval_base_quality_unassessed` | Always present: what RNA reconstruction cannot establish |
| `reverse_complement_support_only`, `mixed_read_orientations` | Orientation of supporting reads relative to the sequencer's input (not RNA strand) |
| `reconstruction_limit:…`, `nondefault_branch_thresholds` | A search limit was reached, or pruning thresholds were changed |

These are prompts for review, not automatic rejections. IDs are full SHA-256
hashes of canonical JSON: candidate IDs hash the event, reference, sequence,
protein and stop status; segment and fragment IDs are scoped to the sample and
source. Hashes are references, not anonymization.

## Compare with predicted proteins

```sh
isovar sv-rna --bam rna.bam --input event.json --output candidates.json \
  --predictions predictions.json
```

`predictions.json` lists proteins predicted from DNA, for example by Varcode:

```json
{"event_id": "GENE1--GENE2", "reference_name": "GRCh38",
 "source": {"annotator": "varcode", "version": "9.4.2", "ensembl_release": 95},
 "predictions": [{"prediction_id": "model-1", "amino_acids": "MAK...", "complete": true}]}
```

`event_id` and `reference_name` must match the event. Each prediction needs a
unique `prediction_id` and `amino_acids`, which may be `null` for an unresolved
prediction. Other fields are kept.

The output's `prediction_comparison` (schema
`isovar.sv_rna_prediction_comparison.v2`) compares every event-linked
annotated-frame translation and exploratory ORF with each prediction. Statuses
distinguish identical sequences from exact fragments and from no match.
A shorter complete ORF is not called a fragment of a longer prediction.
Junction-only support is never reported as whole-protein support.

## Read evidence you can filter on

`record_evidence` holds, for every cited original record, its MAPQ, whether it
has base qualities, lengths and selected native tags. `original_records` holds
the SAM text, and `alignment_metadata` the RG and PG header lines. Absent tags
and MAPQ 255 are `null` in `mapping_quality`; the raw value is kept.

| Field | Meaning |
|---|---|
| `mapping_quality`, `mapping_quality_raw`, `sam_flag` | Alignment placement confidence (null for 255), the original value, and SAM flags |
| `base_qualities_available`, `original_base_qualities_tag_present` | Whether QUAL is stored, and whether an `OQ` tag exists (never substituted for QUAL) |
| `qs`, `dx`, `pi`, `sp` | ONT mean read Q, duplex status (1 duplex, 0 simplex, −1 simplex parent of a duplex), split-read parent and signal offset |
| `NM`, `mg` | Edit distance and gap-compressed identity (%) |
| `AS`, `NH`, `HI`, `nM` | Alignment score, number of alignments and index, STAR mismatches |
| `ms`, `s1`, `s2`, `dv`, `de`, `rl`, `tp` | minimap2 scores, divergence, repetitive seed length and alignment type |
| `rm` | pbmm2 trimmed overlapping query matches |
| `rq`, `np`, `ec` | PacBio predicted accuracy, passes and effective coverage |
| `ic`, `is`, `im` | Iso-Seq consensus inputs and associated reads (not molecule counts) |
| `CB`, `UB`, `XM`, `RG`, `rc`, `ff` | Cell barcode, UMI, Iso-Seq XM, read group, Iso-Seq real-cell flag, CCS failure bits |
| `CR`, `CY`, `UR`, `UY`, `RX`, `QX`, `MI`, `PG` | Raw barcode/UMI bases and qualities, molecular barcode, molecule ID, producing program |

These are native, producer-specific fields, not interchangeable quality scales.
For example, CCS writes `rq=-1` for unpolished reads. To apply your own policy,
join each witness's `observations[id].records` to `record_evidence`, filter, and
recount distinct fragments. Aggregate counts do not survive arbitrary filtering.
`isovar.read_metadata.record_evidence(read)` extracts the same fields from any
pysam record. To filter reads before reconstruction, pass a predicate:

```python
from isovar import ReadCollector

# Keep uniquely aligned reads; keep reads without an NH tag too.
collector = ReadCollector(read_filter=lambda read: not read.has_tag("NH") or read.get_tag("NH") == 1)
```

The result records that a custom filter was used; describe it in
`event_provenance` so the run can be reproduced.

Cell/UMI labels and ONT signal lineage are reported on each junction
(`direct_cell_umi_support`, `direct_read_lineage`) and each ORF's full-interval
support. They are described in [cell/UMI evidence](cell-umi-evidence.md) and
[ONT read lineage](ont-read-lineage.md). Their `complete_label_count` is null
when any supporting segment lacks a resolved label or library.

## Regression data

Synthetic tests cover both strands, spliced and exact breakpoints, homology,
SA-only claims, distant mates, secondary placements, read-group collisions,
missing QUAL, deep noisy coverage, forks and limits. Original reads cover:

- **ATP5MG–KMT2A** (Sid, STAR): reproduces the validated `MAQFVRNLVEKTPALVNG*`
  translation and its noncoding ambiguity. With soft clips, a read-through TruSeq
  adapter makes a second path until the adapter profile trims it. In every Sid
  long-read product this join lies between annotated splice sites of adjacent
  genes on one strand, so RNA alone cannot tell it from read-through.
- **BCR–ABL1** (K562 Iso-Seq): without the SA partners no placed junction is
  invented; with soft clips, breakpoint-clip paths include the validated junction peptides.
- **TPST1–CRCP, FOXO3, PARD3B** (Sid ONT split reads): the breakpoint junction is
  recovered, including PARD3B's 12-nt insertion; no coding frame crosses it.
- **Sid T1 long reads** ([fixtures](../tests/data/fusions/long-read/README.md)):
  PacBio TPST1–CRCP reads with 8 nt of homology (20 fragments), and noisy ONT FOXO3
  reads (12 fragments, 8 cell/UMI labels with unknown library scope).

Out of scope: whole-genome clip realignment, consensus calling, complex
multi-breakend events, automatic adapter-kit inference
([#309](https://github.com/openvax/isovar/issues/309)), initiation inference,
germline/co-somatic attribution ([#297](https://github.com/openvax/isovar/issues/297))
and Vaxrank ranking. Reconciling these paths with Varcode's SV hypotheses is
[#305](https://github.com/openvax/isovar/issues/305).

## Reproduce the osteosarc comparison

From a checkout with Varcode ≥9.4.2 and Ensembl 95 installed:

```sh
python -m examples.osteosarc_sv_validation output-directory
```

This reconstructs TPST1–CRCP, PARD3B–CDKN2B-AS1/CDKN2B and FOXO3–STRADA/CCDC47,
the three strongest named rearrangement leads in the audited catalogue, from
selected original PacBio/ONT records in both orientations, and compares them with
Varcode's predictions from the T1 DNA calls. It writes one reconstruction and
comparison JSON per event and orientation plus `protein_comparison.csv`. To use a
full indexed T1 alignment instead, pass `--bam PATH_OR_URL` (and `--source URL`
for a local copy). Gene-pair labels do not assert a coding fusion.

The fixture builders live with the tests: `tests/data/fusions/build_long_read.py`,
`build_osteosarc.py`, `build_three_fusions.py`, `audit_three_fusions.py`, and
`validation-events.json` with its generator `build_validation_events.py`.

Sources: [SAM](https://samtools.github.io/hts-specs/SAMv1.pdf) and
[SAM tags](https://samtools.github.io/hts-specs/SAMtags.pdf),
[Dorado](https://software-docs.nanoporetech.com/dorado/latest/basecaller/sam_spec/),
[CCS](https://ccs.how/faq/reads-bam.html), [Iso-Seq](https://isoseq.how/isoseq-tags.html),
[minimap2](https://github.com/lh3/minimap2/blob/master/minimap2.1),
[pbmm2](https://github.com/PacificBiosciences/pbmm2),
[fusion reconstruction assessment](https://pmc.ncbi.nlm.nih.gov/articles/PMC6802306/),
[NCBI translation tables](https://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi),
[Kozak 1986](https://pubmed.ncbi.nlm.nih.gov/3943125/).
