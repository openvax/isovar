# Technical-sequence evidence (#302)

Use the existing adapter/tail matcher and immutable ReadSequenceView throughout
read collection, local assembly and SV reconstruction. Do not infer a preparation
from a platform name, or call the verified Sid adapter-containing read a false
fusion. Preserve original sequence, qualities, CIGAR, tags and lineage.

Implement in reviewable stages, each with a version bump and release gates:

1. Typed, explicit library/preparation provenance and feature roles; bounded reuse
   of identical segment annotations across variants and alternative placements.
   Preserve per-RG/mate selection and profile hashes. Add profile FASTA loading
   through the same schema, not a platform subclass hierarchy or an external
   trimming command. Add source-backed example profiles with their limitations.
2. Error-tolerant internal adapter/linker annotation with explicit uncertain
   boundaries and genomic/technical alternatives. Use the existing Edlib matcher
   primitives; internal candidates are not silently excised or concatenated.
   Derived evidence intervals must not span technical joins; retain usable RNA
   on each side with the same parent identity. Feed this policy through the shared
   collector/SV path, rather than adding a special fusion-only classifier.
3. Supply explicit competing genomic interval evidence through the annotation API
   and the aligned-record context. Tied technical/genomic explanations remain
   ambiguous and untrimmed. Keep expensive whole-genome competition separate from
   bounded read annotation. Preserve known upstream tail/segmentation/consensus
   tags and missing metadata; a processed transcript need not expose its adapters.

Tests: unchanged real Sid Illumina R2 and two ONT tail fixtures; processed PacBio
missing-QUAL input; simulated internal concatemer positive clearly identified as
simulated; exact and noisy genomic negatives; both orientations and tied starts;
profile/RG/mate/cache isolation; unknown source metadata; CLI/API parity. Benchmark
bounded cached annotation on the retained original subsets. No uncalibrated
probability, kit guess, invented molecule count or translation claim.

Primary sources inspected 2026-10-29: https://skera.how/ (segmented reads share
parents); https://isoseq.how/isoseq-tags.html (ic/is/im and barcode/UMI metadata);
https://software-docs.nanoporetech.com/dorado/latest/basecaller/read_trimming/
(kit-dependent trimming, orientation and UMI tags, processed sequence limits).
Automatic unknown-kit discovery remains the separately scoped #309; raw-signal
polya estimation is outside this sequence-only change.
