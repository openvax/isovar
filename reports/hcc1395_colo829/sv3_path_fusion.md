# COLO829 SV3: fused original RNA evidence and protein hypotheses

The source candidate BAM was replayed with dense paired-template support,
soft-clipped bases, exon-anchored discovery, and no MAPQ threshold.
The ten pinned original mates of the five exact junction spots were
checked separately. The [machine-readable report](sv3_path_fusion.json)
contains the assembled sequence, spot IDs, path intervals and parameters.

Source: [SRR8615617](https://www.ncbi.nlm.nih.gov/sra/SRR8615617).
Alignment and clipping terms follow the
[SAM specification](https://samtools.github.io/hts-specs/SAMv1.pdf).

The retained joined path has **2,178 bases**,
**3 placed split-junction spots** and
**5 exact 43-base junction spots**.
It has **0 exact CNIH3 exon-anchor bases**
and **0 event-linked frame translations**.
The bounded search retained 80 paths and reported
mate_location_unavailable, path_limit.

| Spot | Exact-window end | Placed split | Junction-end match includes soft clip | Opposite mate on joined path | Opposite mate on CNIH3 exon | MAPQs |
|---|---:|---|---|---|---|---|
| `SRR8615617.101713327` | 2 | yes | yes | yes | no | 1,52,60 |
| `SRR8615617.110737497` | 1,2 | no | no | no | no | 6 |
| `SRR8615617.7312123` | 2 | no | yes | yes | no | 52,60 |
| `SRR8615617.7313798` | 1 | yes | yes | yes | no | 6,52,60 |
| `SRR8615617.97504216` | 2 | yes | yes | yes | no | 1,3,19 |

Mate sequence agreement with the joined path is a paired-fragment
link. It does not fill the unsequenced insert or establish a splice
path. Alternative placements and low MAPQ remain placement uncertainty;
the exact 19-base unplaced join is still directly observed.

| Exploratory ORF | Start-codon spot matches | One end start-to-core join | Start-to-43-mer overlap chain | Full ORF sequence tiling | Start/frame |
|---|---:|---:|---|---|---|
| 37 aa | 47 | 1 | `SRR8615617.109202987` → `SRR8615617.7313798` | 2 distinct spots | unassigned |
| 162 aa | 51 | 2 | `SRR8615617.109202987` → `SRR8615617.7313798` | 7 distinct spots | unassigned |

The two amino-acid sequences come from different ATGs and frames
on the **same assembled nucleotide path**. The tilings use exact
30-base-or-longer overlaps among distinct original spots; the JSON
gives each interval. They support local RNA sequence content, not one
observed full-length molecule or translation initiation.
Start-codon spot matches each have at least 30 exact path bases
covering the ATG. A single read spanning an entire ORF is not required.
Neither start is annotated or experimentally observed. Their shared
junction fragments cannot rank the proteins. There is no observed
mature RNA connection to a versioned CNIH3 exon and no established
frame, so neither exploratory ORF is an SV3 proteoform call.

The input is an alignment-based target-screen fixture, not every RNA
molecule or all possible splice paths. Its search hit the path limit.
The RNA and nominated DNA were obtained from different COLO829
cultures/passages. Matched-normal RNA was not evaluated.
The [prior polished Iso-Seq check](colo829_sv3_polished_check.md)
found no exact SV3 junction in representative sequences. Published
[PacBio original COLO829 subreads](https://downloads-ap.pacbcloud.com/public/dataset/Melanoma2019_IsoSeq/subreads/COLO829T/)
were not ingested: the direct PacBio host returned HTTP 403 and the
Asia-Pacific host timed out on 2026-10-24. Original long-read RNA or an independent splice
witness remains needed to resolve [#563](https://github.com/openvax/isovar/issues/563).

Regenerate with `python -m reports.hcc1395_colo829.build_sv3_path_fusion`.
