# Same-timepoint long/short RNA example: PIP5K1A

Two **complete indexed locus queries**, acquired September 16, 2026, from the
public CC0 osteosarc collection. The query is GRCh38
`chr1:151242177-151242180`, surrounding the source `151242178 AG>A` allele.
Original reads, qualities, CIGARs and tags are unchanged. No downsampling,
alignment editing, realignment or allele-based selection. Primary derivatives
exclude flags 256, 1024 and 2048. Full source URLs, acquisition receipts,
BAM/index SHA256 and original/retained SAM-record hashes are in
`corpus/manifest.json` (these public record hashes are not sequence validation).

| Product | Original regional alignments | Primary Isovar ref/alt/other objects | Protein, assembly on/off |
| --- | ---: | --- | --- |
| UCSF T1 tagged ONT | 95 | 55 / 4 / 1 | 46 / 46 aa |
| UCSF T1 Illumina | 182 | 129 / 1 / 0 | none / none |

The ONT protein passes independent frame, translation and reference-plus-edit
checks using the existing offline Ensembl 87 corpus models. No retained
non-focal CIGAR indels were found in the selected protein's reconstructions
(known frame limitation #265). The Illumina product has one alternate object,
below the unchanged two-read sequence-coverage floor. Its absence of output
is **insufficient alternate support**, not evidence that the variant is absent.

Both are catalogued T1, June 2024, but different processed libraries. These are
not matched technical replicates or proven independent molecules. Cell
composition, depth, alignment and quality affect the comparison; it does not
isolate read length or prove that short reads cannot reconstruct this allele.
Do not pool their counts or infer a malignant-cell assignment from barcodes.
Source definitions: [osteosarc data catalogue](https://osteosarc.com/data/),
[variant](https://osteosarc.com/variant/PIP5K1A-chr1-151242178/).

## Reproduce acquisition

From the repository root, using HTTPS-enabled samtools and Python with pysam:

```python
import gzip, json, tempfile
from pathlib import Path
from tests.data.osteosarc.expansion.acquire import acquire_regions
from tests.data.osteosarc.figure_comparisons.rebuild import rebuild

matrix = json.loads(gzip.decompress(Path(
    "tests/data/osteosarc/expansion/audit/matrix.json.gz").read_bytes()))
variants = [v for v in matrix["variants"] if v["gene"] == "PIP5K1A"]
acquisition = tempfile.mkdtemp(prefix="isovar-pip5k1a-acquisition-")
for source in matrix["sources"]:
    if source["source_id"] in {"53f498a544883d51", "c89442609fffd3f1"}:
        receipt = acquire_regions(source, variants, acquisition, assembly="GRCh38")
        assert receipt["status"] == "ok"
rebuild(acquisition, "/tmp/pip5k1a-new-fixtures")  # destination must not exist
```

The acquisition tool downloads indexes and bounded regions, never the whole
remote BAM. Compare retained record hashes and header/allele identity as well
as archive checksums (BAM compression can vary by tool version).
Tests and figure generation use only the pinned offline derivatives.

## Context, phasing and haplotype extension

`corpus/context-evidence.json.gz` is separately pinned by
`context-manifest.json`. It retains original SAM records for both ZNF436 RNA
products, both CD109 RNA products, and five MAP2 DNA/RNA products, along with
the 35-product MAP2 audit summary and source hashes. Counts refer to products,
not independent cohorts; tagged/deduplicated ONT copies are never added.

- ZNF436: three alternate CB/UMI groups per T1 library, 3 ONT versus 8 Illumina
  read objects; 49 versus 30 aa with assembly. Even the optimistic union of
  short-read flanks is only 113 nt, below the 147 nt needed for 49 aa. The
  matched timepoint does not control library/cell composition.
- CD109: 21 ONT templates carry both alternate alleles, 207 spliced nucleotides
  apart. The short-read product has 22 focal-alt and 19 linked-alt templates,
  but no read/template or CB/UMI observes both. Three original ONT reads exactly
  support the displayed 222-nt two-edit interval, Q20 at every base. Its local
  translation is conditional on the annotated frame, not a whole-protein
  Isovar reconstruction; one 25-mer cannot include both changed codons.
- MAP2: the 22-nt listed deletion, 28-nt deletion alone, and 28-nt deletion
  plus C>A/T>G substitutions are separate explicit hypotheses. Exact anchored
  compound sequence has 24/16 T1/T2 DNA, 4/2 bulk RNA and one deduplicated T2
  ONT template at Q20. No exact isolated 22-nt sequence appears in the audited
  35 products. This is not proof of absence at arbitrary sensitivity.

`build_context.py` packages the read-only diagnostic audit and original
records; it requires `--diagnostics`, `--reference-models` (the original
Ensembl 87 model archive), and `--output`. `verify_context()` independently
recounts the pinned MAP2/CD109 observations, checks the ZNF436 molecule groups
and translations, and checks the 222-nt CD109 witnesses before plotting.
The remaining 30 MAP2 products retain audit summaries/digests, not raw SAMs.

Generate the entire gallery with
`python -m examples.osteosarc_context_figures --output-dir figures/osteosarc`
(install `isovar[plot]`). Original assembly examples come first, followed by
context/haplotype comparisons and explicit fusion RNA windows. The output is
UTC-stamped, with individual white 600-dpi PNG/SVG panels, per-example vector
PDFs, evidence JSON, a page index, and `isovar-all-figures.pdf`.
