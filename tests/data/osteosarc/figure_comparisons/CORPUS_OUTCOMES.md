# Outcome review after the RNA-quality and insertion-boundary fixes

Scope is the current 49-case original-RNA corpus, 16 explicitly labelled stress
cases and the nine-candidate extended footprint set. These are cases/products,
not independent people or molecules, and not every variant on osteosarc.com.
All original BAM records remain unchanged. An unresolved protein is a valid
outcome when the necessary allele, phase or coding evidence is unavailable.

## Original 49 cases

- 40 return a protein agreeing with the nominated edit's reference prediction
  over the recovered context. This includes the newly retained missing-QUAL
  PacBio MT-ND5 case; mitochondrial translation uses the mitochondrial code.
- Two return independently validated RNA proteins different from the single-edit
  prediction: NTF3 (source-linked compound context) and native-GRCh37 NR2F2
  (the RNA-observed nonfocal deletion). Neither result establishes germline or
  somatic origin for the additional change.
- ABCF2, MAP2, MYO15B and MYO9A each have one alternate read object/template
  after collection, below the default two-observation reconstruction floor.
  Overlapping mates are not independent support. No protein is fabricated.
- ADGRF5, the selected GRCh38 NR2F2 product and one MT-ND5 product have no
  callable alternate support. Other products can differ; this is not evidence
  that the variant is biologically absent.

The 1.18.1 insertion-boundary fix leaves all 49 recorded results unchanged.
Every returned protein still passes independent sequence/frame/interval checks.

## Sixteen separate stress cases

Eight ACSL6/EPPK1 product/locus queries have no original aligned records in the
queried regions. Five cases have callable reference/other RNA but no alternate
support. The missing-QUAL PacBio NTF3 compound case now retains one alternate
observation, still below the reconstruction floor. Two ONT NTF3 compound
cases reconstruct the expected compound protein. The only 1.18.1 changes are
count-only corrections in two KTN1 products; all stress proteins are unchanged.

These are separately acquired, coordinate-validated stress queries, not
vaccine-membership claims. Empty queried regions remain distinct from a failed
download, an unmapped allele, or lack of a coding reference.

## Nine extended candidates and SV context

GTF3C5, RNF213, GLIS3 and KTN1 have exact-allele RNA support in at least one
queried product, with all reconstructed alternatives independently translated.
Sample/technology/assembly and missing-QUAL status remain explicit. Alternate
support can be insufficient in an individual product. Lower-ranked GTF3C5
ONT/PacBio discrepancies remain hypotheses, not confirmed biological isoforms.

DLG5 has a sequence-resolved DNA junction, but no mutation-linked RNA/CDS
in the bounded queries. AFF3/KEAP1 overlap normal intronic splice footprints
that do not identify the DNA allele. GABBR1--SLC29A1 and OTUD7A--FMN1 have
observed supplementary RNA paths, without a justified coding reconstruction
from the retained windows. These results are informative, not zero-valued
protein predictions. DNA-nominated WIPF2/TMEM63B clips are retained but do not
yet establish partner assignments. See [the full report](EXTENDED_RNA.md).

The [matched clipping audit](SOFT_CLIPS.md) covers both policies and both
assembly modes; it provides no basis for enabling all soft clips by default.
The next scientific step for unresolved SVs is evidence-linked, splice-aware
partner realignment/consensus and CDS justification, not threshold relaxation
or reference-filled sequence.
