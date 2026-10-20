"""Replay a published canine tumor DNA call against original tumor RNA reads."""

import argparse
import gzip
from pathlib import Path
import shutil
from tempfile import TemporaryDirectory

from pyensembl import Genome
from varcode import load_vcf

from isovar import ProteinSequenceCreator, export_protein_hypotheses, run_isovar, write_protein_hypotheses


DATA = Path(__file__).resolve().parents[1] / "tests" / "data" / "wlos0001_tp53"


def run_example(scratch_dir=None):
    """Return every protein hypothesis from the WLOS_0001 TP53 interval."""
    with TemporaryDirectory(dir=scratch_dir) as directory:
        directory = Path(directory)
        fasta = directory / "canfam3-chr5-tp53-masked.fa"
        with gzip.open(DATA / (fasta.name + ".gz"), "rb") as source, fasta.open("wb") as dest:
            shutil.copyfileobj(source, dest)
        genome = Genome(
            reference_name="CanFam3.1",
            annotation_name="Ensembl",
            annotation_version=98,
            gtf_path_or_url=str(DATA / "tp53-ensembl-98.gtf"),
            transcript_fasta_paths_or_urls=[str(DATA / "tp53-transcripts-ensembl-98.fa")],
            protein_fasta_paths_or_urls=[str(DATA / "tp53-proteins-ensembl-98.fa")],
            genome_fasta_path_or_url=str(fasta),
            cache_directory_path=str(directory / "cache"),
        )
        genome.index()
        variants = load_vcf(str(DATA / "wlos0001-tp53.vcf"), genome=genome)
        results = run_isovar(
            variants, str(DATA / "wlos0001-tp53.bam"),
            protein_sequence_creator=ProteinSequenceCreator(max_protein_sequences_per_variant=0),
        )
        return export_protein_hypotheses(results, sample_id="WLOS_0001", source="SRR8691780")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", default="wlos0001-tp53-protein-hypotheses.json")
    args = parser.parse_args()
    result = run_example()
    write_protein_hypotheses(result, args.output)
    event = result["events"][0]
    print("RNA reads at TP53 5:32563056: ref=%d alt=%d other=%d" % tuple(
        event["allele_support"][allele]["reads"] for allele in ("ref", "alt", "other")))
    print("Protein hypotheses: %d" % len(event["protein_hypotheses"]))
    print("Wrote %s" % args.output)


if __name__ == "__main__":
    main()
