"""Pin Ensembl 87 GRCh38 CNIH3 transcripts for the COLO829 SV3 audit."""

from dataclasses import asdict
import gzip
import json
from pathlib import Path

from pyensembl import EnsemblRelease

from isovar.fusion import FusionReference


ROOT = Path(__file__).resolve().parents[2]
OUTPUT = ROOT / "tests/data/colo829_whole_run/ensembl87_cnih3_models.json.gz"


def build_models(genome=None):
    """Keep CNIH3 cDNA models and only validated complete coding sequences."""
    genome = EnsemblRelease(87) if genome is None else genome
    genes = [gene for gene in genome.genes_by_name("CNIH3") if gene.contig == "1"]
    if len(genes) != 1:
        raise ValueError("Expected one GRCh38 CNIH3 gene")
    gene, = genes
    models, unavailable, skipped = [], {}, {}
    for transcript in sorted(gene.transcripts, key=lambda t: t.id):
        sequence = transcript.sequence.upper()
        if set(sequence) - set("ACGT"):
            skipped[transcript.id] = "ambiguous_cdna_bases"
            continue
        exons = tuple(sorted((exon.start - 1, exon.end) for exon in transcript.exons))
        complete = (transcript.contains_start_codon and transcript.contains_stop_codon
                    and transcript.start_codon_complete and transcript.stop_codon_complete)
        start = transcript.first_start_codon_spliced_offset if complete else None
        end = transcript.last_stop_codon_spliced_offset + 1 if complete else None
        fields = (transcript.id, "GRCh38", "Ensembl 87", "chr1", gene.strand,
                  exons, sequence)
        try:
            reference = FusionReference(*fields, start, end)
        except ValueError:
            unavailable[transcript.id] = "CDS_fails_complete_frame_validation"
            reference = FusionReference(*fields)
        models.append(asdict(reference))
    return dict(schema="isovar.colo829_sv3_reference_models.v1", release=87,
                assembly="GRCh38", gene_id=gene.id, models=models,
                cds_unavailable=unavailable, skipped_transcripts=skipped)


def main():
    data = json.dumps(build_models(), sort_keys=True, separators=(",", ":")).encode()
    with OUTPUT.open("wb") as raw, gzip.GzipFile(fileobj=raw, mode="wb", filename="", mtime=0) as compressed:
        compressed.write(data)
    print(OUTPUT)


if __name__ == "__main__":
    main()
