# Dense exon-anchored SV discovery (#559)

## Problem

The bounded SV search can seed at a clipped end of a versioned reference exon
and bridge through exact observed sequence. Dense discovery currently samples
representatives before this rule runs. Under a record cap, a distant exon clip
can be omitted, and the API rejects the combined options.

## Change

Use the same exon placement, terminal-clip orientation, and minimum clipped
sequence rules in dense representative priority as bounded seeding. Keep the
candidate search and sequence-only extension identical after selection. Count
support from the complete eligible input, preserving original fragment
identity, supplementary placements, raw qualities, and explicit resource
limits. A clip is a candidate seed, never proof of an SV or coding frame.

## Verification

Create a reference-anchored intronic-breakpoint fusion whose only seed is an
exon clip. Include many ordinary reads, an ordinary splice competitor, and
duplicate supplementary payloads. Under a tight discovery cap, dense and
bounded mode must recover the same observed RNA path; dense direct support
must count the distinct originals. Confirm the existing source-derived HCC
fusion replay and the standard lint, test, and CI gates.

The broader #555 work follows in a separate PR: additional sequence-bridge
contexts, source-derived SV3 exon/ORF audit, and regenerated HCC/COLO table.
