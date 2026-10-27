# COLO829 DDX17: why ALT RNA did not yield a protein

The original result was a **coding-reference eligibility failure**, not a negative
test for cryptic splicing. `ProteinSequenceCreator` returns before RNA assembly
when `reference_contexts_for_variant` is empty. DDX17 chr22:38489856 GA>G has
zero such contexts in Ensembl 87. [Isovar #583](https://github.com/openvax/isovar/issues/583)
tracks reusable intronic RNA assessment beyond this report investigation.

## Results

- Original atlas: 221 ALT, 325 REF and 56 OTHER fragments. These groups are
  not mutually exclusive: 11 names occur in multiple allele-observation groups.
  The site-selected genomic BAM retains 679 records for 603 fragment names;
  75 records contain soft clipping and none contains a splice CIGAR operation.
- Recovered all 603 original pairs, including their missing mates. Their original
  SRA placements comprise 870 mapped ends and 336 unplaced ends; no mapped end
  has a splice CIGAR operation. Original placements use GRCh37, not GRCh38.
- Expanded to **27,211 regional pairs plus 98 additional selected pairs = 27,309**
  original pairs, all from SRR8615617. STAR 2.7.11b aligns these to a 50-kb GRCh38
  reference without supplied junction annotations. It reports 25,988 locally
  unique pairs, 433 multiply placed pairs and 888 unaligned pairs. Local uniqueness
  is not genome-wide uniqueness. All placements, including alternatives, are retained.
- The complete 3,940-base intron between coding exons has at least **115 locally
  unique fragments at every base**. This supports intron-containing RNA; tiling
  does not establish one full-length, phased transcript or distinguish retained
  introns from incompletely processed RNA.
- **449 fragments** support the ordinary splice across chr22:38488116-38492055
  (one-based skipped interval). It removes the mutation-bearing intron.
- NMD-model junctions are present, including 12 fragments across
  chr22:38489459-38492055. There are **10 fragments** at the alternative
  intronic acceptor, chr22:38490437-38492055, and two at
  chr22:38491155-38492055. All three have canonical minus-strand splice motifs.
  The 10-fragment junction has 11–48 aligned bases on its shorter anchor.
- Of 221 original ALT fragments, 214 are locally unique, five multiply placed
  and two unaligned. **No locally unique ALT-fragment alignment is spliced.**
  One multiply placed ALT fragment (`SRR8615617.75487984`) has a secondary
  `9M2683N92M` placement across chr22:38486940-38489622. Its 9-base anchor,
  noncanonical motif and absence of annotated exon endpoints leave it ambiguous.
  Nearby noncanonical 10-base-anchor placements also remain in `results.json`;
  they are not silently discarded or promoted to exon-linked evidence.
- Reference coding-frame continuation into the retained intron predicts
  `GLYR*`, stopping at chr22:38492042-38492044. The 10-fragment alternative
  path predicts `EALISAADHSIHGRDKSTLAGSSSR*`, stopping at
  chr22:38490360-38490362. Both complete coding models give these same local
  predictions. **These are reference-path predictions, not inferred ALT proteins.**
  Both stops precede the deletion in minus-strand transcription/translation order.

Ensembl 87 models: ENST00000396821 and ENST00000403230 are complete coding
transcripts; ENST00000216019 is retained-intron and lacks a complete CDS;
ENST00000475004 is an incomplete NMD model. All four full transcript nucleotide sequences
are reproduced exactly by concatenating their exons from the pinned genomic slice.

This audit supports intronic RNA and alternative splicing, but does not establish
a mutation-phased coding path or that the deletion caused a splice change.
Cryptic initiation, other RNA paths and matched-normal comparisons remain unresolved.
The distinction between junction evidence and mutation-bearing junction evidence
follows the [MiSplice primary study](https://www.nature.com/articles/s41467-020-19307-6).

## Sources and selection

`manifest.json` pins every original input, reference and placement by SHA-256.
The only RNA library is [SRR8615617](https://www.ncbi.nlm.nih.gov/sra/SRR8615617),
COLO829 CCLE paired cDNA RNA, distinct from the benchmark-DNA and TAGET cultures.
The original atlas fixture is `tests/data/colo829_whole_run/`.

With NCBI SRA Toolkit 3.4.1, original site-selected fragment numbers were recovered
from the atlas BAM query names and passed to:

```sh
vdb-dump SRR8615617 -T SEQUENCE -R "$selected_spots" -f fastq1
vdb-dump SRR8615617 -T SEQUENCE -R "$selected_spots" \
  -C PRIMARY_ALIGNMENT_ID,READ_LEN -f json -I
vdb-dump SRR8615617 -T PRIMARY_ALIGNMENT -R "$selected_alignment_ids" \
  -C REF_SEQ_ID,REF_POS,REF_ORIENTATION,CIGAR_LONG,SEQ_SPOT_ID,SEQ_READ_ID,MAPQ -f json -I
```

Regional acquisition uses the actual reference index, not a whole-run FASTQ scan:

```sh
vdb-dump SRR8615617 -T REFERENCE --idx-range i_name
vdb-dump SRR8615617 -T REFERENCE -R 573731-573736 \
  -C NAME,SEQ_ID,SEQ_LEN,PRIMARY_ALIGNMENT_IDS -f json -I
```

Only rows **573732 and 573733** were used for the expanded selection. On
[CM000684.1, GRCh37 chr22](https://www.ncbi.nlm.nih.gov/nuccore/CM000684.1),
these represent alignment starts in zero-based `[38880000,38890000)`.
Their 44,668 primary alignment IDs were queried for the metadata preserved in
`regional-original-alignments.json.gz`; distinct `SEQ_SPOT_ID` values yield
`regional-spots.txt`. All reads of those spots, including unplaced mates, were
retrieved with `vdb-dump -T SEQUENCE -R ... -f fastq1`.
Source-alignment selection can miss fragments placed elsewhere or wholly unmapped.
No claim of whole-library completeness is made. Original coordinates are not
used as GRCh38 alignments; the recovered sequence is realigned.

The local reference is GRCh38 chr22:38470001-38520000 from the 1000 Genomes
full-analysis-plus-decoy/HLA FASTA linked in the manifest. STAR's BAM coordinates
are **relative to this slice**, whose contig name embeds the bounds. They must
not be opened against a whole-genome reference without coordinate translation.
`results.json` and the PDF use GRCh38 chromosome coordinates throughout.

## Reproduce

Recompute every reported count, junction, coverage bin and reference-frame prediction:

```sh
python -m reports.hcc1395_colo829.build_ddx17_audit
python -m reports.hcc1395_colo829.build_evidence_pdfs
```

Replay the alignment from the pinned original FASTQs with the official
[STAR 2.7.11b release](https://github.com/alexdobin/STAR/releases/tag/2.7.11b):

```sh
python -m reports.hcc1395_colo829.replay_ddx17_alignment \
  --star /path/to/STAR --work-dir /tmp/ddx17-replay
```

No reference sequence fills unobserved mutant protein in this report. The current
patch improves the report and pins this investigation; it does not implement the
general intronic assembly/reconciliation engine requested in #583.
