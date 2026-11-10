# RNA origin heuristic

Implement the user's requested best-guess NuMT assessment as one public,
category-independent RNA origin API and an `isovar rna-origin` command. Reuse
the original-library index, segment observation builder, path matching and
scoped evidence identities. No separate BAM parser, assembler or genetic-code
resolver. Source labels are explicitly supplied hypotheses, not gene-name or
contig-name inference.

The panel contains nucleotide paths, per-base positions, a nonempty target
interval (allele or anchored junction), and nuclear/mitochondrial origin labels.
Every hypothesis must carry the same target sequence up to reverse complement;
the focal allele must not bias the context comparison. Indel targets include
placed anchor bases. Both origins must be represented. Panels are explicit,
bounded hypotheses: no whole-genome or patient-specific NUMT discovery claim.

For each original segment, identify target-bearing observations by exact target
sequence and coordinate placement. Compare their flanks, anchored at the target,
against every hypothesis with edlib; use the whole eligible observed context,
with free distal reference overhangs only. Require sufficient reference context
for every competitor. Mask low-quality context bases; missing qualities cannot
establish an origin. Initial heuristics: Q20, 30 reliable context bases, maximum
5% context edits for the best fit, and at least three edits of separation between
origin groups. Record all thresholds and scores, and expose them as policy.
These choices are conservative engineering defaults, not calibrated likelihoods.

Retain all original-placement alternatives; they cannot vote independently.
Conflicting origins within a segment or fragment remain conflicting. Uncertain
alternatives cannot be resolved by majority vote. Combine only target-bearing
segments within the same read-group/QNAME fragment. An unobserved/disjoint mate,
SA declaration without an actual record, or overlap assembly across independent
templates cannot establish allele-to-origin phase. Return per-fragment evidence,
distinct scoped segment/fragment sets, panel identity and eligibility provenance.
Different source libraries remain separate. Sample-level mixed evidence remains
mixed. Origin assessment does not select a translation table or remove candidates.

Public API and JSON CLI tests must cover nuclear and mitochondrial context,
reverse orientation, ties, focal-allele neutrality, short panels, low/missing
quality, indel/junction targets, repeated alignments, competing placements, mate
conflicts, read-group/library identity, and offline JSON output. Preserve existing
RNA feature behavior when sharing its coordinate-offset helper. Full lint/tests,
CI, feature-branch PR, version bump and clean-main PyPI deployment are required.

Primary evidence: Wei et al. 2022, https://pmc.ncbi.nlm.nih.gov/articles/PMC9630118/
(nuclear/mtDNA junction linkage); Maude et al. 2019,
https://doi.org/10.3389/fcell.2019.00201 (allele-dependent cross-mapping);
Laricchia et al. 2022, https://pmc.ncbi.nlm.nih.gov/articles/PMC8896463/
(non-reference NUMTs and copy-number-dependent false positives). Native-DNA
methylation is useful independent evidence but is unavailable from ordinary RNA
or cDNA; it is outside this RNA heuristic.
