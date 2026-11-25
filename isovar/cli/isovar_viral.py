"""Recover observed viral RNA sequences from combined-reference alignments."""

from pathlib import Path

from ..default_parameters import GDC_REFERENCE
from ..logging import configure_cli_logging, get_logger
from .help import HelpArgumentParser
from .input_args import add_genome_fasta_arg, add_sample_id_arg, add_scratch_dir_arg
from .output_args import add_force_arg, add_log_level_arg, add_output_dir_arg
from .rna_args import add_rna_path_args, rna_path_context
from .validation import input_error

logger = get_logger(__name__)


def make_parser(prog="isovar viral-rna"):
    parser = HelpArgumentParser(prog=prog, description=__doc__)
    add_rna_path_args(parser)
    add_genome_fasta_arg(parser, required=True, help="Local FASTA used for alignment")
    add_sample_id_arg(parser, required=True)
    add_output_dir_arg(
        parser,
        required=True,
        help="New directory for viral.json, sequences.fa and predicted proteins.fa",
    )
    add_force_arg(parser, help="Allow writing into an existing --output-dir")
    parser.add_argument("--contig", action="append", help="Optional viral contig subset; repeat to select several")
    parser.add_argument("--alignment-reference", choices=[GDC_REFERENCE, "custom"], default=GDC_REFERENCE,
                        help="Checked reference profile (default: %(default)s); custom requires --contig")
    parser.add_argument("--max-discovery-reads", type=int, default=10000,
                        help="Read segments for discovery per virus; support still scans the full input "
                             "(default: %(default)s)")
    parser.add_argument("--max-paths", type=int, default=256, help="Assembled paths per viral contig (default: %(default)s)")
    parser.add_argument("--min-overlap", type=int, default=30, help="Minimum exact overlap for assembly (default: %(default)s)")
    add_scratch_dir_arg(parser, help="Temporary disk-index directory")
    parser.add_argument("--reference-cache", help="PyEnsembl genomic-FASTA cache directory")
    add_log_level_arg(parser)
    return parser


def run(args=None, prog=None):
    parser = make_parser(prog or "isovar viral-rna")
    options = parser.parse_args(args)
    configure_cli_logging(options.log_level)
    try:
        from ..file_identity import file_sha256 as digest
        from ..json_files import publish_json
        from ..viral_coding import write_viral_proteins
        from ..viral_sequences import reconstruct_viral_sequences, write_viral_fasta

        output = Path(options.output_dir).expanduser()
        if output.exists() and not options.force:
            raise ValueError("Output directory already exists: " + str(output))
        fasta = Path(options.genome_fasta).expanduser().resolve()
        bam_path = Path(options.bam).expanduser().resolve()
        reference_name = "GRCh38" if options.alignment_reference == GDC_REFERENCE else None
        with rna_path_context(options, reference_name=reference_name) as (bam, collector, library):
            result = reconstruct_viral_sequences(
                bam, genome_fasta=fasta, sample_id=options.sample_id, source=str(bam_path),
                contigs=options.contig,
                reference_profile=None if options.alignment_reference == "custom" else GDC_REFERENCE,
                max_discovery_reads=options.max_discovery_reads, max_paths=options.max_paths,
                min_overlap=options.min_overlap, scratch_dir=options.scratch_dir,
                reference_cache=options.reference_cache, read_collector=collector, library_index=library)
        result["input"]["sha256"] = digest(bam_path)
        output.mkdir(parents=True, exist_ok=options.force)
        publish_json(output / "viral.json", result)
        write_viral_fasta(result, output / "sequences.fa")
        write_viral_proteins(result, output / "proteins.fa")
        logger.info("Wrote viral.json, sequences.fa and proteins.fa to %s", output)
    except (OSError, ValueError, KeyError, TypeError) as error:
        input_error(parser, error)
