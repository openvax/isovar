# Adapter and poly-A/T inference (#302)

## First-PR design

Add a shared sequence-end annotator, explicit versioned adapter profiles and an
immutable derived sequence view. Annotation and trimming are independent;
existing defaults and original BAM records remain unchanged.

- Match explicit adapters with substitutions/insertions/deletions and terminal
  partial matches. Kit identity is supplied, never guessed from an instrument.
  Profiles can select adapters by mate and can be assigned per read group.
- Infer A/T-rich terminal tails, including interrupted runs. They are candidates,
  not proof of polyadenylation, exact original tail length or transcript completeness.
- Record intervals in the original available SAM query orientation, with a
  reversible map for retained coordinates and original sequencing orientation.
  Hard-clipped bases remain unavailable. No new coordinate scalar type.
- Optional adapter/poly-A trimming returns one contiguous retained interval;
  never stitch across an internal technical join. For aligned reads only
  terminal CIGAR `S` bases are eligible. Preserve all aligned bases and CIGAR `I`
  evidence, including terminal insertions. Do not mutate BAM sequence, qualities,
  CIGAR, MD, SA or read identity.
- Derive alleles from the original alignment, then slice sequence, qualities and
  reference positions together and rebase allele offsets once. Keep original
  coordinates for supplementary-path phasing. Carry optional source-view
  provenance through compact storage, mate merging and allele conversion.
- Expose the same opt-in policy to CLI/API. Unknown quality remains unknown.
  Unknown profiles still permit tail annotation, not guessed adapter removal.

Verify partial/noisy adapters, both orientations, unknown quality, genuine
genomic homopolymers, terminal insertions, hard/soft clips, indels/splices,
supplementary paths, compact/public parity, mate merging and unchanged defaults.
Include unchanged Sid adapter/tail records as regression fixtures.

This is the safe terminal-sequence foundation, not all of #302: automatic kit
identification, raw-signal tail estimation, a comprehensive kit catalogue,
internal concatemer splitting and genomic disambiguation of aligned candidate
sequence remain separate work. Trimming candidate tails is explicitly opt-in.

## Python

```python
from isovar import Adapter, ReadEndProfile, ReadCollector, infer_read_ends

profile = ReadEndProfile(
    name="explicit-library", version="1",
    adapters=(Adapter("R2", "AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT", mate=2),))
collector = ReadCollector(
    use_soft_clipped_bases=True, read_end_profile=profile,
    trim_adapters=True, trim_poly_a=True)
view = collector.read_sequence_view(pysam_read)
print(view.annotations, view.sequence)
# view.original_sequence / original_qualities retain the input unchanged.
# view.original_interval(a, b) maps retained offsets to original SAM SEQ.
# view.sequenced_interval(a, b) also accounts for strand and hard clips.

# Unaligned sequences are supported without a BAM or reference genome:
view = infer_read_ends("GCTACGTCGCTG" + "A" * 26, trim_poly_a=True)
assert view.sequence == "GCTACGTCGCTG"
```

Use `infer_read_ends=True` on `ReadCollector` to annotate without trimming or
selecting a known kit. A profile also enables annotation. For mixed libraries,
pass `read_end_profile={"read-group-ID": profile}`; unlisted groups remain
unknown. Adapter matching is orientation-aware and supports an optional
`both_orientations` profile setting for cDNA libraries. Unknown QUAL is preserved.

Each collected observation's `source_read_views` retains its source alignment
identity and view, including profile name, version and configuration SHA-256.
These survive compact collection and allele conversion.
Merged mates retain separate source views, **not** an invented one-to-one mapping
from the consensus to one original read. Source views precede allele-level N
cleanup and are not coordinates of a merged/assembled sequence. Existing CSV
schemas remain unchanged; inspect annotations through the Python API.

## CLI

All RNA commands accept `--infer-read-ends`, `--read-end-profile FILE.json`,
`--trim-adapters` and `--trim-poly-a`. Trimming defaults to off. To reuse clips
after end trimming, also enable the existing `--use-soft-clipped-bases` option.
Without that option, ordinary soft-clip exclusion still takes precedence.

Example profile (provide the actual library's sequences; do not guess its kit):

```json
{
  "name": "explicit-library",
  "version": "1",
  "adapters": [
    {"name": "R2", "sequence": "AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT", "mate": 2}
  ]
}
```

For mixed read groups, use `{"read_groups": {"RG-ID": PROFILE}}` instead.
Profiles allow `end` (`left`/`right`), `min_overlap` and `max_error_rate` per
adapter, and `end_window`, `min_poly_a_length`, `max_poly_a_error_rate` and
`both_orientations` per library. Thresholds are heuristics, not calibrated
probabilities. Short low-complexity matches, equally scoring conflicting cut
boundaries and window-saturated tail boundaries are not resolved by trimming.
Literal adapter sequence uses A/C/G/T: variable barcode/UMI stretches should not
be passed as informative adapter bases. There is no implicit manufacturer
catalogue. Existing upstream `pt` values, including failure states, are retained
separately and are not used to invent missing sequence or locate a trim boundary.

## Release behavior

Isovar 1.19.0 adds opt-in inference and trimming; default read sequences and
CLI/API policies remain unchanged. No aligned sequence or CIGAR is edited.
Edlib supplies the error-tolerant matcher. Tail classification remains a
sequence heuristic, and an inferred tail does not assign an SV or a transcript.

## Scientific references

- [SAM specification](https://samtools.github.io/hts-specs/SAMv1.pdf): query/CIGAR
  orientation and soft/hard clipping.
- [Cutadapt matching](https://cutadapt.readthedocs.io/en/stable/guide.html#adapter-search-parameters):
  partial matches and overlap-relative error tolerance.
- [Dorado poly-A estimation](https://software-docs.nanoporetech.com/dorado/latest/basecaller/polya_estimation/):
  A/T orientation and estimated lengths distinct from the basecalled sequence.
- [Iso-Seq workflow](https://isoseq.how/clustering/cli-workflow.html): primers and
  tails may have been removed upstream; absence in SEQ is not biological absence.
