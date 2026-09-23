# Three-event original RNA fixtures

The builder and audit are checked in **beside the data**, and the regression
suite is `tests/test_sv_rna_three_fusions.py`. These fixtures preserve the
September 22–23, 2026 investigation of PARD3B/CDKN2B, GABBR1/SLC29A1 and
OTUD7A/FMN1. They support candidate short ORFs, not observed protein expression.

## Rebuild

From the repository root, with Ensembl 115 annotation/cDNA/ncRNA cached:

```sh
python -m tests.data.fusions.build_three_fusions \
  --local-corpus /path/to/dataset \
  --output /tmp/three-fusions
```

The local corpus must contain `source_inventory.json` and indexed
`alignments/<source_id>/reads.bam`. The manifest records each source URL,
local BAM hash, source-header hash, explicit query regions, selection rule,
annotation URLs and fixture hashes. Original inventory entries are retained in
the hash-pinned `source-inventory.json.gz`. Only the four
specified original tagged ONT / mapped PacBio products are used; deduplicated
alternatives are never pooled with their predecessors.

For fresh bounded acquisition, install `isovar[data]`, explicitly create an
osteosarc snapshot containing those assets, then use:

```sh
python -m tests.data.fusions.build_three_fusions \
  --snapshot SNAPSHOT_NAME --cache /path/to/osteosarc-cache \
  --output /tmp/three-fusions
```

This delegates indexed range acquisition to `isovar.sid_data.extract_regions`
and retains its receipts. It does not create a snapshot or download whole
alignment files. Local-corpus and snapshot acquisition receipts differ;
compare the selected SAM checksums and reference models when comparing modes.
If source objects or regional corpus coverage differ, selections may differ:
regeneration is explicit, never an automatic golden-file update during tests.

Reproduce the interpretation, offline from the checked-in fixtures:

```sh
python -m tests.data.fusions.audit_three_fusions \
  --output /tmp/three-fusion-audit.json
./test.sh -q tests/test_sv_rna_three_fusions.py
```

## Construction and verification

1. Event coordinates are 0-based interbase retained donor/acceptor boundaries,
   derived from the [osteosarc fusion catalogue](https://osteosarc.com/fusions/).
   Both ends are queried within +/-2,000 bp. All local original SAM records for
   an RG/QNAME/mate segment observed at both ends are retained; identical SAM
   lines returned by overlapping queries are represented once. This selection
   does not use sequence, quality, CB/UB labels, ORFs or a preferred orientation.
   Records outside the queried windows are not fetched by this recipe.
2. Every SAM line is preserved byte-for-byte in its original textual form and
   keyed by SHA-256. The compact header retains decoding contigs, relevant read
   groups and referenced program ancestry via `minimal_header`; empty fixtures
   retain the queried contigs. This is not the complete producer header and must
   not be used to infer missing library/producer metadata. The source header is
   separately hashed. No read sequence, quality, CIGAR or tag is synthesized.
3. Full Ensembl 115 transcript models for the listed overlapping/nearby genes
   are pinned, including the neighboring MYMX locus and unnamed CDKN2B-region
   transcript. Exons are 0-based half-open; only complete annotated transcripts
   have CDS bounds. The builder skips transcripts without an A/C/G/T sequence.
   These models constrain interpretation; they never fill missing RNA bases.
4. Reconstruction is run in both orientations, with `assemble=False`, default
   read filtering, `min_orf_amino_acids=8`, `max_records=50000`, `max_paths=200`
   and `max_extension_segments=500`. Reversing an event swaps its ends and flips
   both strands. Candidate starts are ATGs, not an exhaustive non-AUG search.
5. Only complete candidate ORFs crossing an actual `breakpoint_junction` enter
   the audit. Full-interval witness identities are unioned across paths. Their
   nucleotides are independently translated by the existing test oracle using
   [NCBI standard code 1](https://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi#SG1).
   Each witness is checked against the **original forward query sequence**,
   separately reporting native versus reverse-complement matches, absent or
   ambiguous primary records, minimum full-interval quality, terminal poly(A),
   and CB/UB pairs. Missing PacBio QUAL is never interpreted as Q0 or Q10.
6. Terminal poly(A) and original query orientation support a proposed RNA
   direction; they are not proof of transcript polarity, initiation or protein
   expression. Internal priming and library orientation remain limitations.
   In particular, `ts` is not independent confirmation when minimap2 used
   `-uf` ([upstream manual](https://lh3.github.io/minimap2/minimap2.html)).
   CB/UB counts are observed labels within one source/RG, not established
   independent molecules or cross-library abundance.

The tests check original-record round trips, exact candidate peptide sequences,
strand-separated witness counts, absence of a claimed annotated junction
translation, missing qualities and retained repeat alternatives. They perform
no network access and do not regenerate or overwrite fixtures.

## Findings reproduced by the tests

| Catalogue label | Candidate in the supported orientation | Exact native full-ORF witnesses |
| --- | --- | --- |
| GABBR1::SLC29A1 | `MKRLVSSSRAWWRMPVIPAPTEAEAGESLESGRRRLQ` (37 aa) | ONT T1: 1; T2: 3; all Q>=10, all terminal A>=10 |
| OTUD7A::FMN1 | `MAAKAGTSLEARSLRPAWPTC` (21 aa), reverse of the catalogue orientation | ONT T2: 9, 7 CB/UB labels; all Q>=10; 8 terminal A>=10 |
| PARD3B::CDKN2B | `MNISNIHISTQKKKKKKSRF` (20 aa), reverse of the catalogue orientation; competing repeat lengths remain | ONT T1: 4; T3: 2 (all Q>=10 and A>=10); PacBio T1: 2 with missing QUAL |

GABBR1 joins to intergenic sequence upstream of SLC29A1; the short candidate
stops before the annotated SLC29A1 gene. The second T1 junction witness contains
the reverse complement of this ORF in its native query and is reported
separately. The six catalogue/local junction reads therefore do not imply six
native 37-aa ORF witnesses. There is one additional T3 junction path, without
this exact complete candidate.

The OTUD7A label is misleading for the observed predominant RNA direction:
FMN1 intronic minus-strand sequence joins OTUD7A antisense plus-strand sequence.
Longer 31–66-aa ORFs obtained in the catalogue orientation have only
reverse-complement witnesses here. The earlier selected fixture had eight
paths; the broader +/-2 kb tagged-BAM selection has eleven segments, nine of
which witness the native 21-aa candidate. These counts have different scopes.

PARD3B similarly joins CDKN2B-side minus-strand RNA into PARD3B antisense
sequence. The 20-aa candidate has cross-platform support, but T1 and T2 also
support competing 17-, 19-, 31- and 32-aa ORFs around the A homopolymer. The
fixture deliberately retains those alternatives. T2 does not supply an exact
full-interval witness for the 20-aa candidate in this analysis. A single noisy
T3 read also yields a longer candidate, with full-interval minimum Q3.

No candidate above uses a supported annotated initiation site. None is an
established conventional coding fusion or demonstrated vaccine antigen.
Annotated-CDS continuation and candidate local ATG initiation are different
claims. Empty local fixtures and missing exact candidates are bounded results,
not proof of biological absence in the full libraries.
