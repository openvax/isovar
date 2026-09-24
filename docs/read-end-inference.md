# Adapter and poly-A/T inference

Isovar can annotate, and optionally trim, adapter and poly-A/T sequence at read
ends. Both are opt-in; annotation and trimming are independent, and original BAM
records are never modified.

## Behavior and limits

- **Adapters** are matched with substitutions, insertions, deletions and terminal
  partial matches. The kit is supplied in an explicit, versioned profile, never
  guessed from the instrument. Profiles can select adapters by mate and can be
  assigned per read group. Every equally optimal adapter start is recovered, and
  competing cut boundaries are preserved rather than trimmed.
- **Poly-A/T tails** are A/T-rich terminal candidates, including interrupted runs.
  A/T support is required in the first `min_poly_a_length` bases and in each
  following window of that length, within the configured error rate; the scan
  stops at unsupported sequence, so a distant tail cannot pull the boundary through
  the read body. A tail is not proof of polyadenylation, of the original tail length
  or of transcript completeness.
- **Coordinates** are recorded in the original SAM query orientation, with a
  reversible map to retained coordinates and to sequencing orientation.
  Hard-clipped bases remain unavailable.
- **Trimming** returns one contiguous retained interval and never stitches across
  an internal technical join. For aligned reads only terminal CIGAR `S` bases are
  eligible; aligned bases and CIGAR `I` evidence, including terminal insertions,
  are kept. Sequence, qualities, CIGAR, MD, SA and read identity are not edited.
- **Alleles** come from the original alignment; sequence, qualities and reference
  positions are then sliced together. Original coordinates are kept for
  supplementary-path phasing, and source-view provenance survives compact storage,
  mate merging and allele conversion. Unknown quality remains unknown.

Not handled: automatic kit identification ([#309](https://github.com/openvax/isovar/issues/309)),
raw-signal tail estimation, a kit catalogue, internal concatemer splitting and
genomic disambiguation of aligned candidate sequence
([#302](https://github.com/openvax/isovar/issues/302)). Edlib supplies the
error-tolerant matcher.

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

## Selection and cost

The ordinary allele collector fetches the requested genomic locus (plus its
one-base anchors), applies read/overlap filters and derives allele coordinates
before optional end inference. It does not preprocess the whole BAM. Retained
reference-supporting and alternate-supporting observations are both eligible;
end inference is not a mutation-support classifier. Adapter/tail searches are
bounded by `end_window` (default 200 bases per end). Full original sequence and
quality provenance is still retained. A read overlapping multiple requested
variants can be processed separately for each locus. Calling the standalone
annotator directly makes the caller responsible for selecting reads.

SV discovery needs a different candidate pool: affected gene/exon regions and
both breakpoint partners, not an exact predicted fusion sequence or only reads
covering the literal DNA breakpoint. [`isovar sv-rna`](sv-rna.md) collects that
pool for one nominated event, and `isovar fusion` validates supplied,
sequence-resolved fusion transcripts. Ordinary splicing alone cannot assign an RNA
footprint to a particular DNA mutation.

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

## Scientific references

- [SAM specification](https://samtools.github.io/hts-specs/SAMv1.pdf): query/CIGAR
  orientation and soft/hard clipping.
- [Cutadapt matching](https://cutadapt.readthedocs.io/en/stable/guide.html#adapter-search-parameters):
  partial matches and overlap-relative error tolerance.
- [Edlib](https://github.com/Martinsos/edlib): infix and anchored-prefix edit-distance
  alignment; tied starts are recovered explicitly before choosing trim boundaries.
- [Dorado poly-A estimation](https://software-docs.nanoporetech.com/dorado/latest/basecaller/polya_estimation/):
  A/T orientation and estimated lengths distinct from the basecalled sequence.
- [Iso-Seq workflow](https://isoseq.how/clustering/cli-workflow.html): primers and
  tails may have been removed upstream; absence in SEQ is not biological absence.
