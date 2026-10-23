# Original-read HCC1395 EIF3K–CYP39A1 fusion fixture

The [public INTEGRATE example](https://sourceforge.net/p/integrate-fusion/wiki/example/)
provides the HCC1395 tumor RNA and tumor/normal DNA subset in
`test-data.tar.gz` (SHA-256
`9a37130e768e424935008f1d0a67050f6e11a242a263ed2f20ba071d94b767cb`).
The [INTEGRATE paper](https://pmc.ncbi.nlm.nih.gov/articles/PMC4691743/)
reports multiple EIF3K–CYP39A1 RNA isoforms. This fixture keeps only the
source reads needed to check one canonical coding path; it is not a sample-wide
RNA or DNA library.

- `original_unmapped.bam` contains all 185 original fusion-enriched Illumina
  records, unchanged from the archive (SHA-256
  `73475e017be92c0ee602f2ed71be64ab1e2822cb99ad43c35db175cf3ee4707b`).
- `coding_tiling.bam` contains 19 original aligned Illumina records selected
  from the archive's `accepted_hits.bam`. Together with split reads they tile
  the complete coding sequence by exact sequence overlaps of at least 20 nt.
- `remapped_grch37_lite.bam` and its index contain the 185 original Illumina
  query sequences realigned to the archive's tailored GRCh37-lite reference
  with minimap2 2.26 (`-ax sr -k11 -w5 -Y`); original query names are preserved.
- `models.json.gz` pins 18 GRCh37-lite exon/cDNA models from the archive.
  Complete CDS offsets from Ensembl 87 are transferred only when the entire
  transcript sequence exactly matches the archive's GRCh37 cDNA (SHA-256
  `6f2b085eaf3f9029a4e750b44c2eb3933b22ae09dc92d9df3d7c6e9d071a2a2d`).
- `ont_candidate_sequences.fasta.gz` contains 21 original query sequences
  from HCC1395 ONT direct RNA [SRR23456709](https://www.ncbi.nlm.nih.gov/sra/SRR23456709).
  They were selected because local GRCh38 minimap2 placements touch both
  genes. Six contain the exact 20+20-base fusion junction sequence. These
  local placements are a candidate screen, not a whole-genome uniqueness test
  (SHA-256 `876639fe516064afb3ef4ad135f0851de4a13849bdc6cdeb2fcde72244e65bf3`).

To regenerate the Illumina files from an extracted archive, use
`python reports/hcc1395_colo829/build_hcc1395_eif3k_fixture.py --source-dir
<archive>/test-data --output-dir <new-fixture-dir>`. Add
`--ont-bam <locally-realigned-SRR23456709.bam>` to regenerate the ONT query
subset. The generator uses the existing local Ensembl 87 cache; tests and the
report use only the pinned models and require no live annotation service.

From the repository root, generate the reviewable sequence and evidence with
`python reports/hcc1395_colo829/build_hcc1395_eif3k_fusion.py --output-dir
<report-dir> --pdf`. The committed example is in
`reports/hcc1395_colo829/eif3k_cyp39a1/`. The 531-amino-acid sequence is a
predicted coding product of this RNA path; neither these reads nor Isovar
measure translated protein expression.
