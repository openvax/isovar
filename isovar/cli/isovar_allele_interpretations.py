"""
Compare competing interpretations of one locus with the RNA reads: which
representation of an SNV, multi-base change or indel do the reads carry?
"""

import argparse

from ..allele_interpretations import interpretations_from_dict, reconcile_allele_interpretations
from ..default_parameters import INTERPRETATION_FLANK, MIN_INTERPRETATION_FRAGMENTS
from ..json_files import read_json, write_json
from ..logging import configure_cli_logging, get_logger
from .output_args import add_log_level_arg
from .rna_args import add_rna_args, alignment_file_from_args, read_collector_from_args
from .validation import CommandInputError, check_outputs, input_error, non_negative_int, positive_int

logger = get_logger(__name__)


def make_parser(prog="isovar allele-interpretations"):
    parser = argparse.ArgumentParser(prog=prog, description=__doc__)
    parser.add_argument("--input", required=True,
                        help="JSON (or .json.gz) with reference_name and interpretations: ID to a list of "
                             "variants with contig, start, ref and alt")
    parser.add_argument("--output", required=True,
                        help="Candidate by RNA evidence JSON, gzip-compressed if it ends in .gz")
    parser.add_argument("--sample-id", required=True, help="Sample the reads came from")
    parser.add_argument("--source", help="Identity of the read set, which scopes read IDs (default: --bam)")
    parser.add_argument("--genome", help="Annotation for the variants (default: the input's reference_name)")
    parser.add_argument("--flank", type=non_negative_int, default=INTERPRETATION_FLANK,
                        help="Exonic bases added on each side of the locus (default: %(default)s)")
    parser.add_argument("--min-fragments", type=positive_int, default=MIN_INTERPRETATION_FRAGMENTS,
                        help="Fragments needed to call an interpretation supported (default: %(default)s)")
    parser.add_argument("--cell-umi-labels", action="store_true",
                        help="Also count UMIs and cells from CB/UB tags (single-cell data)")
    add_rna_args(parser)
    add_log_level_arg(parser)
    return parser


def run(args=None, prog=None):
    parser = make_parser(prog or "isovar allele-interpretations")
    options = parser.parse_args(args)
    configure_cli_logging(options.log_level)
    try:
        check_outputs(options, "output")
        interpretations = interpretations_from_dict(read_json(options.input), genome=options.genome)
        with alignment_file_from_args(options) as alignment_file:
            result = reconcile_allele_interpretations(
                interpretations, alignment_file, sample_id=options.sample_id,
                source=options.source or options.bam, read_collector=read_collector_from_args(options),
                flank=options.flank, min_fragments=options.min_fragments,
                cell_umi_labels=options.cell_umi_labels)
        write_json(options.output, result)
        logger.info("Wrote %s", options.output)
    except (OSError, ValueError, CommandInputError) as error:
        input_error(parser, error)
