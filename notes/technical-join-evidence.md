# Internal technical joins and competing sequence origin (#302)

Build on the typed profiles/cache, preserving all opt-in defaults. Add an explicit
internal motif placement to profile data (separate from terminal end orientation),
source-backed roles and expected adapter order where supplied. Reuse error-tolerant
Edlib matching with full internal motifs and explicit boundary uncertainty; never
call an unknown kit from platform metadata. Keep simulated positive controls
labelled as simulated. Skera documents ordered adapter FASTA and parent/S-read
relationships; processed RNA without visible linkers is a preservation control.

Represent candidate technical intervals over the unchanged original query.
Keep raw qualities, CIGAR, tags, hard clips and parent identity. A derived view
may expose several usable contiguous intervals; it must never concatenate them
across a technical join. Sequence matching remains separate from placement-aware
origin assessment and expensive whole-genome competition.

Add explicit genomic-competitor evidence. Exact supplied/reference-supported
query intervals, mere alignment and unavailable reference context are different
states. Preserve ties and competing explanations; an adapter-like genomic motif
must not be silently deleted. CIGAR insertion bases remain original evidence,
not automatic technical sequence. MD/reference evidence must be interpreted in
SAM query coordinates and verified against the observed sequence; missing MD
is not a negative genomic result. Poly-A/T candidates remain a separate category.

Use the same derived intervals in all consumers:

- Local allele assembly selects the contiguous usable interval containing the
  allele, preserving original quality footprints and source-read views.
- SV observations split into separate intervals with the SAME original identity,
  adjusted positions/breaks/query extent. Do not add synthetic molecules or permit
  self-chaining across a technical boundary. Dense/raw CIGAR shortcuts must apply
  the same policy, not count a join that the built observation excludes.
- Path assessment receives the same collector policy, including raw FASTQ and
  original alternative placements. Exact sequence/exon-clip witnesses must lie
  wholly in an eligible interval. Store original records and annotations for
  review. Routing can overselect raw templates, but evidence scope must pin the
  profile/policy so incompatible reused evidence raises.

Existing exact_path_matches accepts separate sequence intervals with the same
(fragment, mate) identity; interval_spots requires one contiguous end and overlap
chains forbid reuse of a fragment. Reuse that behavior rather than minting mate
numbers or parent-fragment IDs. Original query offsets need explicit preservation
when intervals are sliced. Audit witnesses and counters should explain excluded
technical or unresolved joins without dropping the underlying reads.

Regression tests: internal positive with two biological sides; no source mutation
or false spanning witness; ordinary-reference adapter-like motif; mismatches,
indels and tied boundaries; both orientations and soft/hard clipping; true CIGAR
insertions; same-parent segments cannot become independent coverage; dense versus
exhaustive and routed versus exhaustive policy parity; raw FASTQ; changed-policy
scope validation; unchanged real Sid terminal fixtures. Run a small independent
match oracle for the internal matching algorithm and original-library annotation
benchmarks without claiming whole-library prevalence or translation.

Primary sources read 2026-10-29:
https://skera.how/adapters.html
https://skera.how/read-segments.html
https://skera.how/output.html
https://isoseq.how/isoseq-tags.html
https://software-docs.nanoporetech.com/dorado/latest/basecaller/read_trimming/

Source caution: Skera's read-name description says 0-based [qStart,qEnd), but a
note on the same page says left-open/right-closed. Preserve qs/qe raw and verify
against implementation/PacBio BAM specification before interpreting them. Do not
silently choose a coordinate convention from the conflicting prose.

Implementation decisions in progress: Adapter now has explicit placement
terminal/internal plus min_internal_flank. The internal matcher recovers optimal
overlapping boundaries, then recursively searches disjoint remaining intervals;
a perfect copy cannot suppress a noisier distant copy. Suboptimal overlapping
matches are not counted as additional joins. Terminal and internal matching share
the full Edlib boundary helper. Existing 68 terminal regressions are being checked.

Next: typed original-query genomic evidence (exact reference claim, alignment only,
reference unavailable), fresh origin assessment after cached sequence matching,
and explicit excluded/usable intervals without changing ReadSequenceView.sequence
semantics. A genomic sequence competitor stays an alternative and prevents an
adapter-only cut. Uncertain internal joins can withhold cross-join inference under
the explicit trimming/exclusion policy while preserving raw bases and usable sides.
Pass CIGAR/MD evidence lazily, only when technical matches exist; missing/malformed
MD must remain an unavailable/reference-claim state, not evidence against a genomic
origin. Preserve gene-reference alternatives and document limitations. Poly-A
partial genomic overlap must not break the existing five-aligned-A plus 26-clipped-A
Sid regression. Sources from original alignment tags remain claims, not new FASTA
measurements. Verify any supplied FASTA context separately.

PacBio's primary BAM specification resolves the Skera prose discrepancy:
https://pacbiofileformats.readthedocs.io/en/13.0/BAM.html explicitly defines qs/qe
and segmented QNAME intervals as 0-based, end exclusive. Preserve raw values and
parent provenance; do not equate parent counts with independent RNA molecules.

Self-review found #620: bidirectional SV extension forgot its original segment
choice between directions. A synthetic regression produced AACGTT by joining
opposite ends from two incompatible observations of one parent; independent
parents are a valid positive control. Extension now carries reciprocal observation
choices across directions and, on conflict, tries the opposite order to retain
AACG and CGTT as alternatives under the same total path budget. 199 focused
reconstruction/routing/internal tests pass after the fix; broader gates pending.

Eligibility is deliberately conservative and opt-in: an explicit internal motif
plus trim_adapters makes unresolved candidate intervals barriers to spanning
inference, without removing their raw bases or claiming technical origin. Exact
supplied reference or CIGAR=/MD claims preserve a genomic competitor as usable;
M/I placement alone is labelled a competitor but is insufficient to prove exact
origin. Terminal trimming preserves its existing aligned/insertion safeguards.
No whole-genome negative result is claimed. Missing, malformed or contradictory
MD evidence is explicit; independent reference verification remains distinguishable.
