"""Reconstruct antigen isoforms, then optionally estimate short-read abundance."""

import argparse
import json
from pathlib import Path

from pyensembl import Genome

from ..antigen_reference import (
    GDC_REFERENCE, _annotation_cache, gdc_genome, isoform_targets_from_genome, reference_file_identity, validate_alignment_reference,
)
from ..antigen_isoforms import reconstruct_antigen_isoforms
from ..audit_bundle import digest, publish_json
from ..isoform_quantification import quantify_antigen_isoforms, write_isoform_fastas
from .rna_args import add_rna_path_args, rna_path_context


def make_parser(prog="isovar isoforms"):
    parser = argparse.ArgumentParser(prog=prog, description=__doc__)
    actions = parser.add_subparsers(dest="action", required=True)
    reconstruct = actions.add_parser("reconstruct", help="Reconstruct explicit target genes from RNA alignments")
    add_rna_path_args(reconstruct)
    reconstruct.add_argument("--gene", action="append", required=True, help="Gene ID/name; repeat for each target")
    reconstruct.add_argument("--sample-id", required=True)
    reconstruct.add_argument("--genome-fasta", required=True, help="Local FASTA used for alignment")
    reconstruct.add_argument("--gtf", required=True, help="Local, versioned transcript/CDS annotation")
    reconstruct.add_argument("--annotation-version", required=True)
    reconstruct.add_argument("--annotation-name", default="GENCODE")
    reconstruct.add_argument("--assembly", default="GRCh38")
    reconstruct.add_argument("--alignment-reference", choices=[GDC_REFERENCE, "custom"], default=GDC_REFERENCE,
                             help="Combined-reference profile (default: %(default)s); checked against the BAM")
    reconstruct.add_argument("--output-dir", required=True, help="New result directory")
    reconstruct.add_argument("--reference-cache", help="Reusable annotation-cache root (default: PyEnsembl cache)")
    reconstruct.add_argument("--max-discovery-reads", type=int, default=10000,
                             help="Reads used to discover sequence paths; support still scans the full input")
    reconstruct.add_argument("--max-paths", type=int, default=256)
    reconstruct.add_argument("--min-overlap", type=int, default=30)
    reconstruct.add_argument("--scratch-dir")
    quantify = actions.add_parser("quantify", help="Run installed Salmon against the complete competing transcriptome")
    quantify.add_argument("--input", required=True, help="Reconstruction isoforms.json")
    quantify.add_argument("--transcriptome", required=True, help="Full matching background transcript FASTA")
    quantify.add_argument("--reads1", required=True, help="Short-read FASTQ, or first mates")
    quantify.add_argument("--reads2", help="Second mates of paired short reads")
    quantify.add_argument("--output-dir", required=True, help="New directory for Salmon inputs, estimates and report")
    quantify.add_argument("--threads", type=int, default=1)
    quantify.add_argument("--salmon", default="salmon", help="Installed Salmon executable")
    return parser


def run(args=None, prog=None):
    parser = make_parser(prog or "isovar isoforms")
    options = parser.parse_args(args)
    try:
        if options.action == "quantify":
            result = json.loads(Path(options.input).read_text())
            quantify_antigen_isoforms(result, transcriptome=options.transcriptome, reads1=options.reads1,
                                     reads2=options.reads2, output_dir=options.output_dir,
                                     salmon=options.salmon, threads=options.threads)
            return
        output = Path(options.output_dir).expanduser()
        output.mkdir(parents=True, exist_ok=False)
        if options.alignment_reference == GDC_REFERENCE:
            if options.assembly != "GRCh38":
                raise ValueError("The GDC profile uses GRCh38; use custom for other assemblies")
            genome = gdc_genome(genome_fasta=options.genome_fasta, gtf=options.gtf,
                                annotation_version=options.annotation_version, annotation_name=options.annotation_name,
                                cache_directory=options.reference_cache)
        else:
            genome = Genome(options.assembly, options.annotation_name, options.annotation_version,
                            gtf_path_or_url=str(Path(options.gtf).expanduser().resolve()),
                            genome_fasta_path_or_url=str(Path(options.genome_fasta).expanduser().resolve()),
                            cache_directory_path=_annotation_cache(options.gtf, directory=options.reference_cache))
        with genome, rna_path_context(options, reference_name=options.assembly) as (bam, collector, library):
            genome.index()
            targets, unavailable = isoform_targets_from_genome(genome, options.gene, contig_names=bam.references)
            identity = validate_alignment_reference(
                bam, genome, profile=None if options.alignment_reference == "custom" else GDC_REFERENCE,
                contigs=[t.reference.contig for t in targets])
            identity["files"] = reference_file_identity(options.genome_fasta, options.gtf)
            result = reconstruct_antigen_isoforms(
                bam, targets=targets, sample_id=options.sample_id, source=str(Path(options.bam).resolve()),
                reference_identity=identity, max_discovery_reads=options.max_discovery_reads,
                max_paths=options.max_paths, min_overlap=options.min_overlap, scratch_dir=options.scratch_dir,
                read_collector=collector, library_index=library)
            result["input"]["sha256"] = digest(options.bam)
            result["target_selection"] = dict(genes=options.gene, unavailable_transcripts=unavailable)
        publish_json(output / "isoforms.json", result)
        write_isoform_fastas(result, output / "transcripts.fa", output / "proteins.fa")
    except (OSError, ValueError, KeyError, TypeError) as error:
        parser.error(str(error))
