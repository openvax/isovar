"""Recover observed viral RNA sequences from combined-reference alignments."""

import argparse
from pathlib import Path

import pysam

from ..antigen_reference import GDC_REFERENCE
from ..audit_bundle import digest, publish_json
from ..viral_sequences import reconstruct_viral_sequences, write_viral_fasta
from ..viral_coding import write_viral_proteins


def make_parser(prog="isovar viral"):
    parser = argparse.ArgumentParser(prog=prog, description=__doc__)
    parser.add_argument("--bam", required=True, help="RNA aligned to the combined human/viral reference")
    parser.add_argument("--genome-fasta", required=True, help="Local FASTA used for alignment")
    parser.add_argument("--sample-id", required=True)
    parser.add_argument("--output-dir", required=True, help="New directory for viral.json, sequences.fa and predicted proteins.fa")
    parser.add_argument("--contig", action="append", help="Optional viral contig subset; repeat to select several")
    parser.add_argument("--alignment-reference", choices=[GDC_REFERENCE, "custom"], default=GDC_REFERENCE,
                        help="Checked reference profile (default: %(default)s); custom requires --contig")
    parser.add_argument("--max-discovery-reads", type=int, default=10000,
                        help="Read segments for discovery per virus; support still scans the full input")
    parser.add_argument("--max-paths", type=int, default=256, help="Assembled paths per viral contig")
    parser.add_argument("--min-overlap", type=int, default=30, help="Minimum exact overlap for assembly")
    parser.add_argument("--scratch-dir", help="Temporary disk-index directory")
    parser.add_argument("--reference-cache", help="PyEnsembl genomic-FASTA cache directory")
    return parser


def run(args=None, prog=None):
    parser = make_parser(prog or "isovar viral")
    options = parser.parse_args(args)
    try:
        output = Path(options.output_dir).expanduser()
        if output.exists():
            raise ValueError("Output directory already exists: " + str(output))
        fasta = Path(options.genome_fasta).expanduser().resolve()
        bam_path = Path(options.bam).expanduser().resolve()
        with pysam.AlignmentFile(str(bam_path), reference_filename=str(fasta)) as bam:
            result = reconstruct_viral_sequences(
                bam, genome_fasta=fasta, sample_id=options.sample_id, source=str(bam_path),
                contigs=options.contig,
                reference_profile=None if options.alignment_reference == "custom" else GDC_REFERENCE,
                max_discovery_reads=options.max_discovery_reads, max_paths=options.max_paths,
                min_overlap=options.min_overlap, scratch_dir=options.scratch_dir,
                reference_cache=options.reference_cache)
        result["input"]["sha256"] = digest(bam_path)
        output.mkdir(parents=True, exist_ok=False)
        publish_json(output / "viral.json", result)
        write_viral_fasta(result, output / "sequences.fa")
        write_viral_proteins(result, output / "proteins.fa")
    except (OSError, ValueError, KeyError) as error:
        parser.error(str(error))
