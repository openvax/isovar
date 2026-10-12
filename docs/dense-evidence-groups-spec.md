# Dense SV evidence grouping (#464)

Reduce repeated computation in deeply covered regions without inferring cell
origin or treating similar reads as the same molecule. Keep every original SAM
record in the disk spool and count candidate support over the complete input.

## Design

1. Separate computational structure from evidence identity. Bin segment bundles
   by their complete actual alignment structure, then subdivide by full observed
   sequence. Names, cell/UMI labels and quality strings do not create new sequence
   alleles. Different placements, coverage, insertions, clipping, supplementary
   pieces and sequence differences remain distinct. No local-window fingerprint
   or approximate consensus may erase a distant allele or bridge missing bases.
2. Automatically group reads in the existing full-input (`dense_support`) mode;
   introduce no new tuning parameter. Select deterministic original representatives of each complete
   sequence/placement group, retaining enough distinct fragments for the existing
   discovery threshold before spending the rest of `max_records` on original
   voters. Keep abundance: representative-only sampling can promote a rare allele
   by flattening major/minor ratios. A bounded geometry cache shares coordinates
   across same-placement reads, including different alleles, but retains each
   read's sequence and qualities. Keep `max_records` as the existing allocation ceiling.
   The bounded regional mode is unchanged. Report selected
   and omitted groups, records and the representative policy explicitly. Discovery
   votes and quality-dependent annotations remain scoped to the chosen originals.
3. Reject impossible support witnesses using CIGAR joins and breakpoint clips
   before constructing per-base observations or cell/lineage metadata. Actual
   multi-piece/alternative alignments receive the existing full validation; SA
   tags alone do not create a join. Reuse calculations only when the complete
   computational inputs agree, and retain original qualities and contributor
   identities when producing evidence.
4. Benchmark grouping and complete reconstruction separately on pinned short-
   and long-read inputs. Compare candidate sequences and complete-input support,
   runtime, peak memory and scratch use. A grouping compression ratio is not a
   demonstrated reconstruction speedup. Preserve failed/unfinished IGKC attempts;
   their separately budgeted recovery remains #470.

## Verification

Regression cases cover abundant copies and rare late events, substitutions and
indels within the same structural bin, different coverage and placements,
quality differences/missing qualities, missing supplementary pieces, mates,
conflicting cell/UMI labels, repeated signal ancestry, ordering and budget limits.
Check CLI and replay parameter round trips. Run lint and the full test suite,
publish a versioned PR, wait for green CI, then merge and deploy under the
repository release workflow.

Alignment semantics follow the [SAM specification](https://samtools.github.io/hts-specs/SAMv1.pdf)
and [SAM tags specification](https://samtools.github.io/hts-specs/SAMtags.pdf).
Grouping is computational reuse, not molecular deduplication, error correction,
cell-type assignment, or proof of somatic origin.
