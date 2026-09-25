"""
Compare competing interpretations of one locus with the RNA reads: which
representation of an SNV, multi-base change or indel do the reads carry?
"""

import argparse
import json
from pathlib import Path

from ..allele_interpretations import interpretations_from_dict, reconcile_allele_interpretations
from ..default_parameters import INTERPRETATION_FLANK, MIN_INTERPRETATION_FRAGMENTS
from ..logging import configure_cli_logging
from .output_args import add_log_level_arg
from .rna_args import add_rna_args, alignment_file_from_args, read_collector_from_args
from .validation import CommandInputError, non_negative_int, positive_int


def make_parser(prog="isovar allele-interpretations"):
    parser = argparse.ArgumentParser(prog=prog, description=__doc__)
    parser.add_argument("--input", required=True,
                        help="JSON with reference_name and interpretations: ID to a list of variants")
    parser.add_argument("--output", required=True, help="Candidate by RNA evidence JSON")
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
        interpretations = interpretations_from_dict(
            json.loads(Path(options.input).expanduser().read_text()), genome=options.genome)
        with alignment_file_from_args(options) as alignment_file:
            result = reconcile_allele_interpretations(
                interpretations, alignment_file, sample_id=options.sample_id,
                source=options.source or options.bam, read_collector=read_collector_from_args(options),
                flank=options.flank, min_fragments=options.min_fragments,
                cell_umi_labels=options.cell_umi_labels)
        Path(options.output).expanduser().write_text(json.dumps(result, indent=2) + "\n")
    except (OSError, ValueError, CommandInputError) as error:
        parser.error(str(error))
