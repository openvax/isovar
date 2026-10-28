# Viral coding annotation validation — 2026-10-15

Isovar 1.57.0 consumes ORFdata 0.29.0's sequence-qualified HPV16 coding
references, with no copied CDS definitions in Isovar. The upstream source and
offline release are documented in [ORFdata PR #99](https://github.com/openvax/orfdata/pull/99).
The nine source translations, joined coding paths and native/GDC placements
are checked upstream. Isovar's DNA-only fixture retains the source identity in
`tests/data/viral/hpv16-reference.json`.

The consumer regression tests cover the observed sample E7 substitution, both
source-defined joined products, partial ends, reverse-oriented RNA, a mismatching
reference and CLI protein export without a GTF. Synthetic nonviral coding
fixtures check insertion/frame changes, altered initiation and premature stops.
Predictions retain all reference frame anchors and distinguish altered paths
from compatible source paths. A complete prediction does not establish a full
transcript or translation in the sample.

## Saved real RNA replay

The bounded original-read regression now lives in
[`tests/data/viral/siha`](../tests/data/viral/siha/README.md). Its 23-fragment
slice retains the same selected coding prediction, with an offline CLI test,
original SAM hashes, reference/catalogue pins and a deterministic extraction
recipe. The unchanged full-input receipt is stored separately; its 71-fragment
counts below are not counts from the bounded fixture.

The existing public SiHa-1 RNA experiment uses GSM4785628 / SRR12650977,
project PRJNA663801. Its one-million-pair prefix and competitive human/viral
alignment are described in [the recovery experiment](https://github.com/openvax/isovar/issues/514).
No reads were downloaded again or modified for this replay.

The command was rerun on its saved `direct-competitive.bam` with the full GDC
GRCh38.d1.vd1 FASTA, selecting HPV16, with discovery limits of 500 segments and
32 paths. The source FASTA SHA256 is
`fbde121e898f00cc064eaafe5acc58940d7fe069ad14d34acd40c5d675e6b69b`.
The report retains actual BAM/reference checksums and ORFdata catalogue receipts.

- The same **26 RNA paths** and **71 target fragments** were retained. Every
  preexisting per-path field except the newly populated coding fields matched
  the earlier nucleotide-only reconstruction exactly, including all sequence,
  placement and support information.
- These paths yielded **33 coding hypotheses**, including partial products and
  source ambiguity. Counts across hypotheses must not be summed as independent
  RNA or protein abundance.
- The 344-nt path `isovar_viral_a22b9eba3846fcf4496e065f` contains an observed
  start-to-stop, **98-aa E7 prediction** with the sample's L28F difference from
  NP_041326.1. The coding path matches the source placement; the predicted
  protein differs from the reference. The sample sequence is not distributed
  as a replacement reference in ORFdata.
- All nine source CDSs were available for annotation, but only products
  overlapping the observed RNA received predictions. Missing coding ends and
  unassigned fragments remained unresolved.

## Read recovery choice

The head-to-head recovery experiment favored direct minimap2 `splice:sr` with
`-K 1m`, preserving its tested alignments while reducing query-batch memory.
A separate k-mer screen did not improve the measured runtime/sensitivity tradeoff.
This release documents that external workflow without adding an alignment or
screening framework. Human-reference competition remains necessary for
interpreting panel hits; read recovery and RNA/protein annotation are separate
steps.
