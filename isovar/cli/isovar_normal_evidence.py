"""Assess matched-normal genotypes and reference confidence from a pinned inventory."""

import argparse
import logging
from pathlib import Path

from ..normal_sources import run_normal_evidence_manifest
from .commands import parser_for_program


parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("--manifest", required=True, type=Path,
                    help="Versioned normal source inventory and pinned literal-allele nominations")
parser.add_argument("--output", required=True, type=Path,
                    help="New JSON or JSON.gz evidence receipt; existing receipts are preserved")


def run(args=None, *, prog=None):
    command = parser_for_program(parser, prog)
    options = command.parse_args(args)
    logging.basicConfig(level=logging.INFO, format="%(message)s")
    try:
        run_normal_evidence_manifest(options.manifest, options.output)
    except (OSError, ValueError, TypeError, KeyError) as error:
        command.error(str(error))
    print(options.output)
