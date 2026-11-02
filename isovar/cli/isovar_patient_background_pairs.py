"""Export paired mutant and patient-background coding windows from original RNA."""

import gzip
import json
from pathlib import Path

from ..logging import configure_cli_logging, get_logger
from ..patient_background import BASELINES, export_patient_background_pairs
from .commands import parser_for_program
from .main_args import make_isovar_arg_parser, run_isovar_from_parsed_args
from .output_args import add_log_level_arg
from .validation import check_outputs, input_error

logger = get_logger(__name__)

parser = make_isovar_arg_parser(description=__doc__)
parser.set_defaults(max_protein_sequences_per_variant=0)
add_log_level_arg(parser)
parser.add_argument('--sample-id', required=True, help='Sample the reads came from')
parser.add_argument('--source', help='Original RNA read-set identity; defaults to --bam')
parser.add_argument('--baseline', choices=BASELINES, default='focal',
                    help='Somatic edits to revert; germline and unexplained edits are retained '
                         '(default: %(default)s)')
parser.add_argument('--output', required=True, help='Paired-window JSON or JSON.gz output')


def run(args=None, *, prog=None):
    command = parser_for_program(parser, prog)
    options = command.parse_args(args)
    configure_cli_logging(options.log_level)
    try:
        output = Path(options.output)
        check_outputs(options, 'output')
        result = export_patient_background_pairs(
            run_isovar_from_parsed_args(options), sample_id=options.sample_id,
            source=options.source or options.bam, baseline=options.baseline)
        with (gzip.open if output.suffix == '.gz' else open)(output, 'wt') as handle:
            json.dump(result, handle, sort_keys=True, indent=2)
            handle.write('\n')
        logger.info('Wrote %s', output)
    except (OSError, ValueError) as error:
        input_error(command, error)
