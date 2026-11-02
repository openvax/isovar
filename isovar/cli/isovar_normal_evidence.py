"""Assess matched-normal genotypes and reference confidence from a pinned inventory."""

import argparse
from pathlib import Path

from .output_args import add_log_level_arg
from ..logging import configure_cli_logging
from ..normal_sources import run_normal_evidence_manifest
from .commands import parser_for_program
from .validation import input_error


parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("--manifest", required=True, type=Path,
                    help="Versioned normal source inventory and pinned literal-allele nominations")
parser.add_argument("--output", required=True, type=Path,
                    help="New JSON or JSON.gz evidence receipt; existing receipts are preserved")

add_log_level_arg(parser)


def run(args=None, *, prog=None):
    command = parser_for_program(parser, prog)
    options = command.parse_args(args)
    configure_cli_logging(options.log_level)
    try:
        run_normal_evidence_manifest(options.manifest, options.output)
    except (OSError, ValueError, TypeError, KeyError) as error:
        input_error(command, error)
    print(options.output)
