"""Reconstruct antigen isoforms, then optionally estimate short-read abundance."""

import argparse
import json
from pathlib import Path

from gtfparse import ParsingError

from .output_args import add_log_level_arg
from ..logging import configure_cli_logging, get_logger
from ..antigen_reference import GDC_REFERENCE
from ..antigen_targets import CTA_GENE_SET, prepare_antigen_reference
from ..antigen_isoforms import reconstruct_antigen_isoforms
from ..audit_bundle import digest, publish_json
from ..isoform_quantification import quantify_antigen_isoforms, write_isoform_fastas
from .rna_args import add_rna_path_args, rna_path_context
from .validation import input_error

logger = get_logger(__name__)


def make_parser(prog="isovar isoforms"):
    parser = argparse.ArgumentParser(prog=prog, description=__doc__)
    actions = parser.add_subparsers(dest="action", required=True)
    reconstruct = actions.add_parser("reconstruct", help="Reconstruct selected genes or a named catalogue from RNA alignments")
    add_rna_path_args(reconstruct)
    selection = reconstruct.add_mutually_exclusive_group(required=True)
    selection.add_argument("--gene", action="append", help="Gene ID/name; repeat for each target")
    selection.add_argument("--gene-set", choices=[CTA_GENE_SET], help="Versioned OncoRef CTA catalogue")
    reconstruct.add_argument("--sample-id", required=True, help="Sample the reads came from")
    reconstruct.add_argument("--genome-fasta", required=True, help="Local FASTA used for alignment")
    reconstruct.add_argument("--gtf", required=True, help="Local, versioned transcript/CDS annotation")
    reconstruct.add_argument("--annotation-version", required=True,
                             help="Release of the --gtf annotation, e.g. 115 for Ensembl or 47 for GENCODE")
    reconstruct.add_argument("--annotation-name", default="GENCODE",
                             help="Source of the --gtf annotation, e.g. Ensembl (default: %(default)s)")
    reconstruct.add_argument("--assembly", default="GRCh38",
                             help="Assembly of the reads, FASTA and GTF (default: %(default)s)")
    reconstruct.add_argument("--alignment-reference", choices=[GDC_REFERENCE, "custom"], default=GDC_REFERENCE,
                             help="Combined-reference profile (default: %(default)s); checked against the BAM")
    reconstruct.add_argument("--output-dir", required=True, help="New result directory")
    reconstruct.add_argument("--reference-cache", help="Reusable annotation-cache root (default: PyEnsembl cache)")
    reconstruct.add_argument("--max-discovery-reads", type=int, default=10000,
                             help="Reads used to discover sequence paths; support still scans the full input "
                                  "(default: %(default)s)")
    reconstruct.add_argument("--max-paths", type=int, default=256,
                             help="Assembled paths kept per discovery scope (default: %(default)s)")
    reconstruct.add_argument("--discovery-partition", choices=['panel', 'gene'], default='panel',
                             help="Apply discovery budgets across the panel or independently per gene "
                                  "(default: %(default)s)")
    reconstruct.add_argument("--min-overlap", type=int, default=30,
                             help="Minimum exact overlap for assembling reads (default: %(default)s)")
    reconstruct.add_argument("--scratch-dir", help="Existing directory for temporary disk indexes")
    quantify = actions.add_parser("quantify", help="Run installed Salmon against the complete competing transcriptome")
    quantify.add_argument("--input", required=True, help="Reconstruction isoforms.json")
    quantify.add_argument("--transcriptome", required=True, help="Full matching background transcript FASTA")
    quantify.add_argument("--reads1", required=True, help="Short-read FASTQ, or first mates")
    quantify.add_argument("--reads2", help="Second mates of paired short reads")
    quantify.add_argument("--output-dir", required=True, help="New directory for Salmon inputs, estimates and report")
    quantify.add_argument("--threads", type=int, default=1, help="Salmon threads (default: %(default)s)")
    quantify.add_argument("--salmon", default="salmon", help="Installed Salmon executable (default: %(default)s)")
    for action in (reconstruct, quantify):
        add_log_level_arg(action)
    return parser


def run(args=None, prog=None):
    parser = make_parser(prog or "isovar isoforms")
    options = parser.parse_args(args)
    configure_cli_logging(options.log_level)
    try:
        if options.action == "quantify":
            result = json.loads(Path(options.input).read_text())
            quantify_antigen_isoforms(result, transcriptome=options.transcriptome, reads1=options.reads1,
                                     reads2=options.reads2, output_dir=options.output_dir,
                                     salmon=options.salmon, threads=options.threads)
            return
        output = Path(options.output_dir).expanduser()
        # Created only after reconstruction succeeds, so a failed run can simply be repeated.
        if output.exists():
            raise ValueError("Output directory already exists: %s" % output)
        with rna_path_context(options, reference_name=options.assembly) as (bam, collector, library):
            prepared = prepare_antigen_reference(
                bam, genes=options.gene, gene_set=options.gene_set,
                genome_fasta=options.genome_fasta, gtf=options.gtf,
                annotation_version=options.annotation_version, annotation_name=options.annotation_name,
                assembly=options.assembly, alignment_reference=options.alignment_reference,
                reference_cache=options.reference_cache)
            result = reconstruct_antigen_isoforms(
                bam, targets=prepared.targets, sample_id=options.sample_id, source=str(Path(options.bam).resolve()),
                reference_identity=prepared.identity, max_discovery_reads=options.max_discovery_reads,
                max_paths=options.max_paths, min_overlap=options.min_overlap, scratch_dir=options.scratch_dir,
                discovery_partition=options.discovery_partition,
                read_collector=collector, library_index=library)
            result["input"]["sha256"] = digest(options.bam)
            result["target_selection"] = prepared.target_selection
        output.mkdir(parents=True, exist_ok=False)
        publish_json(output / "isoforms.json", result)
        write_isoform_fastas(result, output / "transcripts.fa", output / "proteins.fa")
        logger.info("Wrote isoforms.json, transcripts.fa and proteins.fa to %s", output)
    except (OSError, ValueError, KeyError, TypeError, ParsingError) as error:
        input_error(parser, error)
