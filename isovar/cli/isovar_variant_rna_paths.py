"""Reconstruct original ALT/REF RNA paths before assessing supplied CDS models."""

from contextlib import ExitStack
from pathlib import Path

from ..logging import configure_cli_logging, get_logger
from .help import HelpArgumentParser
from .input_args import (
    add_genome_fasta_arg,
    add_input_arg,
    add_library_cache_arg,
    add_scratch_dir_arg,
    add_source_arg,
)
from .output_args import add_force_arg, add_log_level_arg, add_output_file_arg
from .rna_args import add_rna_args, read_collector_from_args
from .validation import input_error

logger = get_logger(__name__)


def make_parser(prog='isovar variant-rna-paths'):
    parser = HelpArgumentParser(prog=prog, description=__doc__)
    add_rna_args(parser)
    parser.set_defaults(min_mapping_quality=0)
    add_input_arg(parser, required=True, help='Explicit literal variant, assembly and versioned models JSON')
    add_source_arg(parser, required=True, help='Original RNA library identity; keep cultures/libraries separate')
    add_output_file_arg(parser, required=True, help='New JSON or JSON.gz output file')
    add_force_arg(parser, help='Overwrite an existing --output file')
    add_library_cache_arg(parser)
    add_scratch_dir_arg(parser, help='Existing directory for temporary disk indexes')
    add_genome_fasta_arg(
        parser,
        dest='reference_fasta',
        aliases=('--reference-fasta',),
        help='Reference FASTA for decoding CRAM input',
    )
    for name, default in [('max-discovery-reads', 10000), ('max-seed-segments', 10000), ('max-paths', 256),
                          ('min-overlap', 30), ('min-anchor-bases', 12), ('flank', 1000)]:
        parser.add_argument('--' + name, type=int, default=default, help='Default: %(default)s')
    add_log_level_arg(parser)
    return parser


def run(args=None, *, prog=None):
    parser = make_parser(prog or 'isovar variant-rna-paths')
    options = parser.parse_args(args)
    configure_cli_logging(options.log_level)
    try:
        import pysam
        from ..contig_resolution import variant_reference_name
        from ..json_files import read_json
        from ..rna_library import prepare_rna_library
        from ..sv_rna_audit import write_json
        from ..variant_rna_paths import reconstruct_variant_rna_paths, variant_rna_path_input_from_dict

        output = Path(options.output)
        if output.exists() and not options.force:
            raise ValueError('Output exists; choose a new result path')
        supplied = variant_rna_path_input_from_dict(read_json(options.input))
        collector = read_collector_from_args(options)
        parameters = {name: getattr(options, name) for name in
                      ('max_discovery_reads', 'max_seed_segments', 'max_paths', 'min_overlap', 'min_anchor_bases', 'flank')}
        with ExitStack() as contexts:
            bam = contexts.enter_context(pysam.AlignmentFile(options.bam, reference_filename=options.reference_fasta,
                                         threads=options.num_rna_decompression_threads))
            library = None
            if options.library_cache:
                prepared = prepare_rna_library(options.bam, options.library_cache, read_collector=collector,
                                               reference_name=variant_reference_name(supplied['variant']))
                library = contexts.enter_context(prepared.open(
                    mitochondrion_by_length=collector.match_mitochondrion_by_length))
            result = reconstruct_variant_rna_paths(bam, **supplied, **parameters, source=options.source,
                read_collector=collector, library_index=library, scratch_dir=options.scratch_dir)
        output.parent.mkdir(parents=True, exist_ok=True)
        write_json(output, result)
        logger.info("Wrote %s", output)
    except (OSError, ValueError, TypeError) as error:
        input_error(parser, error)
