# Preserve reference competitors to candidate tails

Fix #638 through the existing read-origin model, without introducing a second
tail classifier. Apply placement evidence to every annotation. Before an
opt-in tail cut, assess the actual removed query interval after protecting
aligned bases. Any supplied exact or partial placement competitor prevents
that cut; absent context remains unavailable, not evidence of a biological tail.
Sequence matching alone remains cacheable, and current placement evidence is
evaluated afresh. Original sequence, qualities, CIGAR and tags remain unchanged.

Regression coverage: A/T and both ends; exact and partial supplied competitors;
reference breaks; aligned body plus unaligned tail; missing evidence; cached
sequence with changing evidence; callable laziness and invalid bounds. Keep the
API/CLI policy and documentation consistent. Run focused and full tests, lint,
review the final diff, then exact-head CI, merge and clean-master deployment.

This is an evidence-preservation fix, not calibrated biological-tail inference.
ONT documents that cDNA can contain either A or T tails, failed `pt` estimates
are not tail absence, and estimates do not rewrite the basecalled sequence:
https://software-docs.nanoporetech.com/dorado/latest/basecaller/polya_estimation/
PacBio documents upstream primer/tail removal, so processed inputs need not
contain visible tails: https://isoseq.how/clustering/cli-workflow.html
