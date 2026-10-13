# Execute dense reconstruction once per equivalent observation

Follow-up to #464 and #475. Target an order-of-magnitude reduction in total
reconstruction time, including input handling, discovery and complete support.
Measure against released 1.49.0; do not infer speedup from compression ratios.

1. Profile the pinned short-read IGKC and long-read workloads to identify
   repeated calculations and input/metadata costs. Preserve comparable baselines.
2. Replace expensive repeated operations with shared immutable computation and
   compact group membership where the actual inputs are equivalent. Preserve
   full sequence/placement differences, original records and qualities, fragment
   multiplicity, label resolution, lineage and search-limit semantics. No read
   dropping, approximate consensus, biological deduplication or new tuning knobs.
3. Reduce repeated input parsing, temporary storage and candidate-support work.
   Keep memory bounded by the existing search policy and stream background.
   Validate any broader restructuring against complete baseline results.
4. Add regression cases for unequal allele abundance, quality/label differences,
   mates, supplementary pieces and shared observations across candidate paths.
   Benchmark end-to-end short/long reads and capped high-depth IGKC, recording
   runtime, memory, scratch and candidate/support equivalence. Investigate any
   changed scientific result before shipment.
5. Run lint and the complete test suite, bump the version, create the PR, obtain
   green CI, merge and deploy from clean master, and verify the live release.

If profiling contradicts the proposed design, revise this spec before proceeding.
Report achieved speedups and remaining bottlenecks, including any missed target.

## Profile-guided implementation

The 1.49.0 IGKC profile spends 107 seconds in extension, 100 seconds building
observations (including rebuilding support witnesses), 55 seconds finding and
ranking overlapping reads, and 26 seconds spooling SAM text. These instrumented
times overlap and are not end-to-end speedup measurements.

- Rank identical placement groups, retaining the same original segment cap and
  identity tie-breaks. Group per-base extension votes with original membership.
- Share frozen junction summaries and exact interval comparisons across paths;
  expand successful witnesses to their original query intervals and identities.
- Reuse support observations and bound geometry reuse by both entries and bases.
- Index original BAM offsets and normalize repeated discovery signatures into
  an allele table. Preserve all alternate placements, mates and supplementary
  records; convert SAM/CRAM once without dropping original records.
- Skip ancestry record lookup only when the declared platform already excludes
  that producer model. Preserve original tags and unknown ancestry status.

Alignment and BGZF virtual-offset semantics follow the
[SAM specification](https://samtools.github.io/hts-specs/SAMv1.pdf) and
[SAM tags](https://samtools.github.io/hts-specs/SAMtags.pdf).

Further profiling motivates low-compression scratch for converted inputs and
integer array arithmetic for exact block coverage. Preserve overlapping
alternative-block multiplicity, original segment caps and ties.

Replace collation and BAM payload copying with virtual offsets in the original
BAM. Hash only records sharing a template (the only possible exact duplicates),
visit them in physical order, and retain distant mates/pieces via indexed
lookup. SAM/CRAM inputs require one private BAM conversion. Inputs remain
read-only. Private-input diagnostics remain outside published benchmark data.

Move overlap ranking into integer arrays too. Rank groups by the same coverage
and directional extent, then merge original IDs within tied ranks. A tied group
whose minimum ID follows at least the remaining cap's worth of other minima
cannot contribute a selected identity; partial selection bounds that merge.
