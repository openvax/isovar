# Soft-clipped RNA: impact audit

`use_soft_clipped_bases=False` excludes unaligned read ends, not CIGAR
insertions or deletions. Tests cover a 128-nt insertion and deletion with
both clipping policies and both locus-storage representations. Hard-clipped
sequence is absent from that SAM record and cannot be recovered by this flag.
These distinctions follow the [SAM specification](https://samtools.github.io/hts-specs/SAMv1.pdf).

The flag does **not** realign clips, identify a rearrangement partner, nominate
an allele, or infer a fusion CDS. If a large indel is represented only by a
clip, the small-variant collector does not recover its alternate allele merely
by retaining that clip. A separate one-sided insertion-boundary bug was
reproduced and filed as [#296](https://github.com/openvax/isovar/issues/296):
such reads can incorrectly count as reference support with either setting.
This audit does not fix or silently reinterpret that classification.

## Matched on/off experiment

The same four exact indels, eleven RNA products, Ensembl 87 annotation and
Isovar 1.18.0 branch were run with clipping off/on. All other settings were
identical: absent QUAL retained explicitly as unknown, secondary alignments
excluded, default balanced protein ranking, uncapped protein results, and
assembly separately off/on. This is 44 input comparisons / 88 mode-specific
protein comparisons, not 88 independent biological replicates. It does not
estimate the effect on all SVs or indel sizes.

Reference and alternate RG/QNAME template counts were unchanged in all 44
comparisons. One other-allele count changed (GTF3C5 T0 BostonGene: 3 to 4);
retaining clips can change mate merging even without changing focal bases.
The top protein and full ranked sequence list changed in 22/88 comparisons,
covering eleven product/locus combinations in both assembly modes. All changed
top contexts were shorter; no previously absent protein was recovered.
Every returned protein in both runs passed independent translation, frame
and mutation-interval checks. This validates arithmetic/sequence consistency,
not the biological truth of clipped sequence.

Representative assembly-on top lengths:

| Variant / RNA product | Clips off | Clips on |
| --- | ---: | ---: |
| GTF3C5 / T1 Illumina | 49 aa | 29 aa |
| RNF213 / T1 Illumina | 49 aa | 25 aa |
| GLIS3 / T1 Illumina | 46 aa | 25 aa |
| KTN1 / T2 ONT | 29 aa | 25 aa |
| GTF3C5 / T1 PacBio | 49 aa | 49 aa |

GLIS3 explains why more input sequence is not automatically more usable
context. With clips removed, eight templates are sequence-compatible with
the 46-aa candidate. Retaining the clips introduces conflicting ends and
reduces its compatible template support to one; the same 46-aa protein
remains as rank 3. A 25-aa candidate with eight templates becomes rank 1.
One mate pair also no longer merges. No new frame is established: the shared
protein context agrees. This is the interaction of clip conflicts with the
existing support-aware ranking, not evidence for a shorter biological protein.

## Consequence for SV reconstruction

Keep the small-variant default off. Preserve original BAM sequence and clips
for a breakpoint-aware analysis, then use observed supplementary alignments
or realign clipped sequence to DNA-nominated retained partner/splice paths.
Require consistent orientation, query adjacency and allele linkage, test
competing placements, and retain ambiguous assignments. Only then should
those bases contribute to a mutation-assigned reconstruction/CDS. Missing
QUAL can remain unknown; it must not erase otherwise usable sequence or be
replaced by MAPQ-derived base scores.

The existing fusion-window extractor already reads full original query
sequence and validates actual same-segment supplementary mappings; it does
not depend on this small-variant flag. Unmapped clips, paired mates and
alternative placements are not interchangeable with that path evidence.

## Reproduction and retained inputs

The dated gallery retains the original bounded BAM/index slices, acquisition
receipts, both complete JSON audits and a machine-readable comparison.
From a checkout with the same annotation:

```python
import json
from pathlib import Path
from tests.data.osteosarc.figure_comparisons.footprints import audit

inputs = Path("/path/to/retained/rna")
sources = json.loads((inputs / "sources.json").read_text())
for clips in (False, True):
    audit(inputs, Path("/new/audit-clips-%s" % clips), sources=sources,
          all_proteins=True, use_soft_clipped_bases=clips)
```

Use product `counts` for this primary-only experiment. `default_counts` is
the separate production-default comparator and intentionally does not inherit
the experimental clipping setting.
