"""Export paired mutant and patient-background coding windows from original RNA."""

import gzip
import json
from pathlib import Path

from ..logging import configure_cli_logging
from ..patient_background import BASELINES, export_patient_background_pairs
from .commands import parser_for_program
from .main_args import make_isovar_arg_parser, run_isovar_from_parsed_args
from .output_args import add_log_level_arg
from .validation import check_output_path

parser = make_isovar_arg_parser(description=__doc__)
parser.set_defaults(max_protein_sequences_per_variant=0)
add_log_level_arg(parser)
parser.add_argument('--sample-id', required=True)
parser.add_argument('--source', help='Original RNA read-set identity; defaults to --bam')
parser.add_argument('--baseline', choices=BASELINES, default='focal',
                    help='Somatic edits to revert; germline and unexplained edits are retained')
parser.add_argument('--output', required=True, help='Paired-window JSON or JSON.gz output')


def run(args=None, *, prog=None):
    command = parser_for_program(parser, prog)
    options = command.parse_args(args)
    configure_cli_logging(options.log_level)
    try:
        output = Path(options.output)
        inputs = list(options.vcf or []) + [getattr(options, name, None) for name in
                                          ('bam', 'germline_vcf', 'phasing_error_rates')]
        if output.resolve() in {Path(path).resolve() for path in inputs if path}:
            raise ValueError('Output must not overwrite an input file')
        check_output_path(output)
        result = export_patient_background_pairs(
            run_isovar_from_parsed_args(options), sample_id=options.sample_id,
            source=options.source or options.bam, baseline=options.baseline)
        with (gzip.open if output.suffix == '.gz' else open)(output, 'wt') as handle:
            json.dump(result, handle, sort_keys=True, indent=2)
            handle.write('\n')
    except (OSError, ValueError) as error:
        command.error(str(error))
