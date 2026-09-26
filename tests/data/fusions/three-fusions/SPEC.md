Preserve the construction of the three-event September 2026 RNA investigation alongside regression tests.

- Acquire or reuse bounded original tagged ONT T1/T2/T3 and mapped PacBio T1 records for PARD3B/CDKN2B, GABBR1/SLC29A1 and OTUD7A/FMN1.
- Select segments touching both breakpoint windows without selecting a favored ORF; retain original SAM fields and hashes. Pin full Ensembl 115 transcript models. (The construction code was retired in 1.39.3; the reads are openvax-v1 members, and the builder is at tag v1.39.2.)
- Run reconstruction in both nominated orientations. Independently translate each complete candidate and audit native query orientation, full-interval quality and CB/UB labels.
- Assert reproducible original-read examples, including opposite-strand witnesses, missing PacBio qualities, and competing PARD3B repeat lengths. Candidate ATGs remain unproven initiation sites.
- Document the selection rule (re-applied by the tests), bounded-selection limitations and scientific findings. No production inference changes or vaccine changes.
