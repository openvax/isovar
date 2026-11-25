"""Assess matched-normal genotypes and reference confidence from a pinned inventory."""

from pathlib import Path

from .help import HelpArgumentParser
from .input_args import add_manifest_arg
from .output_args import add_force_arg, add_log_level_arg, add_output_file_arg
from ..logging import configure_cli_logging
from .commands import parser_for_program
from .validation import input_error


parser = HelpArgumentParser(description=__doc__)
add_manifest_arg(
    parser,
    required=True,
    type=Path,
    help="Versioned normal source inventory and pinned literal-allele nominations",
)
add_output_file_arg(
    parser,
    required=True,
    type=Path,
    help="New JSON or JSON.gz evidence receipt file; existing receipts are preserved unless --force is set",
)
add_force_arg(parser, help="Overwrite an existing --output receipt file")

add_log_level_arg(parser)


def run(args=None, *, prog=None):
    command = parser_for_program(parser, prog)
    options = command.parse_args(args)
    configure_cli_logging(options.log_level)
    try:
        from ..normal_sources import run_normal_evidence_manifest

        run_normal_evidence_manifest(options.manifest, options.output, force=options.force)
    except (OSError, ValueError, TypeError, KeyError) as error:
        input_error(command, error)
    print(options.output)
