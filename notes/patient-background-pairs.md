# Patient-background peptide pairs (#100)

Add a native-results export and matching `isovar patient-background-pairs` CLI.
Keep one record per event, nucleotide translation and compatible transcript.
Replay attributed transcript/cDNA edits exactly before reverting either the
focal somatic edit or all supplied somatic edits observed on that RNA. Retain
germline and unexplained changes, and retain co-somatic changes under the focal
policy. Report the actual removed and retained edits explicitly.

Compare bounded coding windows from the same observed RNA extent. Re-translate
after nucleotide edits; allow unequal peptide lengths and unobserved stops.
Never extend with annotation sequence beyond the reverted edit. Ambiguous
edit replay or frame boundaries produce explicit unavailable records.

Computed baselines are counterfactuals, not patient-normal observations. Measure
separate exact reference-allele RNA coverage and spanning-original-read support
using the existing coding coverage function. Report phase as direct spanning
RNA, allele-observation span, assembly-only or unavailable; tiling is not a
single molecule. Preserve original evidence identities and incomplete contexts.

Cover SNVs, germline indels, focal frameshifts/complex indels, minus strand,
alternative transcript coordinates, local caps, uncertain phase and no output.
Exercise the API and CLI using the same construction and schema; document the
scope and downstream comparison policy. Version 1.64.0.

Scientific basis: proximal germline and somatic changes alter the peptide
sequence and require explicit phasing rather than amino-acid substitution:
https://pubmed.ncbi.nlm.nih.gov/30510237/ (Hundal et al., Nature Genetics 2019).
VCF phase sets describe supplied genotype phase, not proof that a local RNA
assembly or reference-allele observation came from normal tissue:
https://samtools.github.io/hts-specs/VCFv4.5.pdf.
