"""
Export every protein hypothesis for each variant, with the nucleotide
translations behind it and scoped RNA evidence, as JSON plus a TSV with one
row per translation. All proteins are kept unless
--max-protein-sequences-per-variant says otherwise.
"""

import sys

from ..logging import configure_cli_logging
from ..protein_hypotheses import export_protein_hypotheses, write_protein_hypotheses
from .commands import parser_for_program
from .main_args import make_isovar_arg_parser, run_isovar_from_parsed_args
from .output_args import add_log_level_arg
from .rna_args import alignment_file_from_args
from .validation import CommandInputError, check_output_path

parser = make_isovar_arg_parser(description=__doc__)
parser.set_defaults(max_protein_sequences_per_variant=0)
add_log_level_arg(parser)
_group = parser.add_argument_group("Protein hypothesis export")
_group.add_argument("--sample-id", required=True, help="Sample the reads came from")
_group.add_argument(
    "--source",
    help="Identity of the read set, which scopes read IDs with --sample-id (default: --bam)")
_group.add_argument(
    "--output", default="isovar-protein-hypotheses.json",
    help="JSON output; the TSV is written beside it (default: %(default)s)")


def run(args=None, *, prog=None):
    command = parser_for_program(parser, prog)
    args = command.parse_args(sys.argv[1:] if args is None else args)
    configure_cli_logging(args.log_level)
    try:
        check_output_path(args.output)
        results = run_isovar_from_parsed_args(args)
        with alignment_file_from_args(args) as alignment_file:
            header = alignment_file.header.to_dict()
        export = export_protein_hypotheses(
            results, sample_id=args.sample_id, source=args.source or args.bam, alignment_header=header)
        write_protein_hypotheses(export, args.output)
    except (CommandInputError, ValueError) as error:
        command.error(str(error))
