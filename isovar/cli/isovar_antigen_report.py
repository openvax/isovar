"""Combine source-specific Isovar exports into a visual antigen evidence audit."""

import argparse
import json
from pathlib import Path

from .output_args import add_log_level_arg
from ..logging import configure_cli_logging
from ..antigen_report import compose
from .commands import parser_for_program
from .validation import CommandInputError, check_outputs, input_error


parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("--sample-id", required=True,
                    help="Sample every input must name; inputs for another sample are rejected")
parser.add_argument("--output", required=True, type=Path, help="Output PDF")
parser.add_argument("--summary-json", type=Path,
                    help="Full source-separated summary JSON (default: --output with .json)")
parser.add_argument("--small-variant-tsv", type=Path,
                    help="Per-call TSV from isovar evidence-report")
parser.add_argument("--sv-result", type=Path, action="append", default=[],
                    help="Versioned isovar sv-rna JSON; repeat per nominated event")
parser.add_argument("--fusion-result", type=Path, action="append", default=[],
                    help="Versioned isovar fusion JSON; repeat per supplied transcript")
parser.add_argument("--viral-result", type=Path, action="append", default=[],
                    help="Versioned isovar viral JSON; repeat per RNA panel/run")
parser.add_argument("--isoforms-result", type=Path, action="append", default=[],
                    help="Versioned isovar isoforms JSON; repeat per target panel/run")
add_log_level_arg(parser)


def run(args=None, *, prog=None):
    command = parser_for_program(parser, prog)
    parsed = command.parse_args(args)
    configure_cli_logging(parsed.log_level)
    if parsed.output.suffix.lower() != ".pdf":
        command.error("--output must end in .pdf")
    summary_json = parsed.summary_json or parsed.output.with_suffix(".json")
    try:
        for path in (parsed.output, summary_json):
            path.parent.mkdir(parents=True, exist_ok=True)
        check_outputs(parsed, "output", summary_json)
        inventory = compose(
            parsed.sample_id, small_variants=parsed.small_variant_tsv,
            sv=parsed.sv_result, fusion=parsed.fusion_result,
            viral=parsed.viral_result, isoforms=parsed.isoforms_result)
        from ..antigen_evidence_pdf import write_antigen_report
        write_antigen_report(parsed.output, inventory)
        summary_json.write_text(json.dumps(inventory, indent=2, sort_keys=True) + "\n")
    except (CommandInputError, OSError, ValueError, KeyError, TypeError) as error:
        input_error(command, error)
    print(parsed.output)
    print(summary_json)
