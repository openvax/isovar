# Cohort report scope and plot readability (#623, #624)

Preserve the difference between an assessed zero ncORF overlap and an absent
catalogue. Count per-nomination catalogue assessment from native states, report
partial coverage explicitly, and retain unknown scope for older results. Reuse
one display helper for count tables and a separate factual scope sentence; do
not change scientific reconstruction or inflate the mutant-protein funnel.

For plots, remove default minor log labels and choose a small, readable set of
integer endpoint/intermediate counts. Keep inclusive survival data unchanged.
Place the funnel legend outside its data area. Test narrow/singleton/broad depth
ranges and mixed catalogue scope, then render and visually inspect the original
DDX17/GOLGA5 report. Preserve every input/result checksum when rendering saved
results. Bump patch version, run required checks and CI, then merge/release.
