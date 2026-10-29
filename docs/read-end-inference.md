# Adapter and poly-A trimming

Read ends can contain technical adapters or primers, and biological poly-A/T
tail candidates. These are distinct annotations: a soft clip, a high base quality
or an A-rich run does not establish technical origin. Isovar marks supplied
technical motifs and candidate tails over unchanged reads and can optionally
trim eligible terminal sequence from a derived view.

- Both are **off by default**, and marking (annotation) and trimming are
  separate choices.
- Only terminal soft-clipped bases can be trimmed. Aligned bases and CIGAR
  insertions are always kept.
- The BAM records themselves are never modified.

## Quick start

Describe your library's adapters in a profile. Use the actual sequences from
the kit you used; Isovar does not guess the kit:

```json
{
  "name": "explicit-library",
  "version": "1",
  "adapters": [
    {"name": "R2", "sequence": "AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT", "mate": 2}
  ]
}
```

Then run any RNA command with it:

```sh
isovar run --vcf variants.vcf --bam rna.bam --output results.csv \
    --read-end-profile library.json --trim-adapters --trim-poly-a \
    --use-soft-clipped-bases
```

| Option | Effect |
|---|---|
| `--infer-read-ends` | Annotate candidate adapter and poly-A/T ends without a profile; nothing is trimmed |
| `--read-end-profile FILE.json` | Use this library's adapters; also turns on annotation |
| `--trim-adapters` | Trim matched adapters from terminal soft clips; needs a profile |
| `--trim-poly-a` | Trim candidate poly-A/T tails from terminal soft clips; turns on annotation |

Soft clips are still excluded unless you also pass `--use-soft-clipped-bases`.
With it, whatever remains of a clip after trimming is used.

## The profile

| Setting | Level | Meaning |
|---|---|---|
| `name`, `version` | library | Identify the profile; recorded with every annotated read |
| `adapters[].name`, `sequence` | adapter | The literal technical motif, in A/C/G/T |
| `adapters[].role`, `source` | adapter | `adapter`, `primer` or `linker`, and its sequence source |
| `provenance` | library | Supplied platform/instrument, preparation, kit/revision, library state, processing steps and source receipts |
| `adapters[].mate` | adapter | Match only on mate 1 or 2 |
| `adapters[].end` | adapter | Match only at the `left` or `right` end |
| `adapters[].min_overlap`, `max_error_rate` | adapter | How short and how noisy a match may be |
| `end_window` | library | Bases searched at each read end (default 200) |
| `min_poly_a_length`, `max_poly_a_error_rate` | library | How long and how clean a tail must be |
| `both_orientations` | library | Search both read orientations, for cDNA libraries whose reads can come from either strand |

For libraries mixing read groups, map read group IDs to profiles:
`{"read_groups": {"RG-ID": PROFILE}}`. Read groups not listed get the empty
`unknown` profile, so no adapters are searched in them.

Keep variable barcode and UMI stretches out of adapter sequences, since they are
not informative adapter bases. There is no built-in catalogue of manufacturer
kits. The thresholds are heuristics, not calibrated probabilities.

## Source provenance and FASTA profiles

`ReadLibraryProvenance` records supplied facts, with unknown fields left as null.
`ReadProcessingStep(tool, operation, version, source)` records upstream processing
without inferring that every read was successfully trimmed. Neither `PL` nor an
instrument name selects a kit. Provenance and motif roles enter the profile
fingerprint; existing inline JSON profiles remain supported.

A profile can instead reference a FASTA file relative to the JSON file:

```json
{
  "name": "custom-preparation", "version": "1",
  "provenance": {"preparation": "custom", "sources": ["sample-sheet.tsv"]},
  "adapter_fasta": {
    "path": "primers.fa", "sha256": "SHA256_OF_ORIGINAL_FASTA_BYTES",
    "defaults": {"role": "primer", "end": "right", "mate": 2,
                 "min_overlap": 12, "max_error_rate": 0.1}
  }
}
```

The first FASTA header token becomes the motif name. Duplicate IDs, ambiguous
bases, empty inputs, FASTQ and changed checksums are errors. The example's
thresholds are explicit configuration, not a calibration for a named kit. All
entries share `defaults`; use inline adapters for differing mate/end policies.
The same declaration works inside a per-read-group profile map. The Python API
is `read_end_profile_from_fasta(path, sha256=digest, profile=fields,
adapter_defaults=settings)`. A linker role here labels terminal matches; internal
concatemer handling remains a separate step in #302.

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

`ReadCollector(infer_read_ends=True)` annotates without trimming or a profile.
For mixed libraries, pass `read_end_profile={"read-group-ID": profile}`.

## How ends are recognized

- **Adapters** are matched allowing substitutions, insertions, deletions and a
  partial match at the read end. Adapter matching is orientation-aware. Every
  equally good adapter start is found. If they imply different cut points,
  the read is left untrimmed.
- **Poly-A/T tails** are A/T-rich runs at a read end, interruptions allowed. The
  first `min_poly_a_length` bases, and each following window of that length,
  must be A/T within the error rate. The scan stops at the first window that
  fails, so a distant tail cannot pull the cut point through the read body.
  A tail does not prove polyadenylation, the original tail length or that the
  transcript is complete.
- **Trimming** keeps one contiguous interval of each read and never stitches
  across an internal technical join.
- **Alleles** are read from the original alignment. Sequence, qualities and
  reference positions are then sliced together, and missing qualities stay
  missing. The read's sequence, qualities, CIGAR, MD, SA tags and identity are
  not edited.
- **Coordinates** are recorded in the original SAM orientation, with a
  reversible map to the retained sequence and to sequencing orientation.
  Hard-clipped bases stay unavailable.

Some cases are deliberately left untrimmed: short low-complexity matches,
conflicting cut points with equal scores, and tails that fill the whole search
window. Upstream `pt` tags, including failure states, are kept separately and
are not used to invent sequence or place a cut.

## Which reads are examined

Only the reads collected for each requested variant are examined, after the
usual read and overlap filters; Isovar does not preprocess the whole BAM.
Reads supporting the reference and the alternate allele are treated alike;
end inference does not classify mutation support. Searches look at most
`end_window` bases from each end. A collector reuses exact sequence/profile/mate/orientation matches across variants
and placements. Its cache holds at most 1,024 entries and one million sequence
bases. Original identity and lineage remain attached separately; cached matches
never add molecules. Each view applies the current alignment's qualities, CIGAR
bounds, hard clips and upstream tags afresh. Direct `infer_read_ends` calls can
supply a `ReadEndInferenceCache` explicitly, with configurable bounds.

## Where the annotations go

Each collected observation's `source_read_views` records the source alignment
and its view, with the profile name, version and configuration SHA-256. These
survive compact collection and conversion to alleles. Merged mates keep one view
per mate, not an invented mapping from the merged sequence to one original read.
Views are recorded before allele-level N cleanup and do not give coordinates in
a merged or assembled sequence. The CSV outputs do not include annotations; read
them through the Python API.

## Not handled

- Identifying an unknown kit automatically ([#309](https://github.com/openvax/isovar/issues/309))
- Estimating tails from raw signal
- Splitting internal concatemers
- Deciding whether an aligned candidate sequence is genomic
  ([#302](https://github.com/openvax/isovar/issues/302))

## References

- [SAM specification](https://samtools.github.io/hts-specs/SAMv1.pdf): query and
  CIGAR orientation, soft and hard clipping.
- [Cutadapt matching](https://cutadapt.readthedocs.io/en/stable/guide.html#adapter-search-parameters):
  partial matches and overlap-relative error tolerance.
- [Edlib](https://github.com/Martinsos/edlib), which supplies the error-tolerant
  matcher: infix and anchored-prefix edit-distance alignment. Tied starts are
  recovered explicitly before trim points are chosen.
- [Dorado poly-A estimation](https://software-docs.nanoporetech.com/dorado/latest/basecaller/polya_estimation/):
  A/T orientation, and estimated lengths distinct from the basecalled sequence.
- [Iso-Seq workflow](https://isoseq.how/clustering/cli-workflow.html): primers and
  tails may have been removed upstream, so their absence in SEQ is not
  biological absence.

## Cohort and SV audit manifests

Use the same inline profile definition under `read_collection.read_end_profile`,
or reference a profile JSON through the manifest's checksum-verified `files`:

```json
{"read_collection": {"read_end_profile": {"file_id": "library-profile"},
                     "trim_adapters": true}}
```

The `library-profile` file entry uses the normal `path` and `sha256` fields.
Any FASTA path inside that JSON is relative to the JSON file and has its own
required digest. Both profile loading and FASTA verification occur before resume
reuses results. The cohort task's input identity includes the profile file.

Cohort tasks already select one source. SV audit sources may supply their own
`read_collection` settings, overriding manifest defaults for that source only.
Use separate source profiles when libraries have different preparations, even if
they reuse a read-group name. An omitted override inherits the explicit manifest
default; no metadata inference selects another kit. Unknown read groups within a
profile map remain unknown.

`ReadCollector.from_dict(settings, files=verified_paths)` is the shared Python
entry point for these serialized settings. The normal constructor continues to
accept typed profiles and mappings of RG to typed profile. Neither loader mutates
the caller's settings.
