"""Combine source-specific Isovar exports into a visual antigen evidence audit."""

import argparse
import json
from pathlib import Path

from ..antigen_report import compose
from .commands import parser_for_program


parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("--sample-id", required=True)
parser.add_argument("--output", required=True, type=Path, help="Output PDF")
parser.add_argument("--summary-json", type=Path, help="Full source-separated summary JSON")
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


def run(args=None, *, prog=None):
    parsed = parser_for_program(parser, prog).parse_args(args)
    if parsed.output.suffix.lower() != ".pdf":
        parser.error("--output must end in .pdf")
    try:
        inventory = compose(
            parsed.sample_id, small_variants=parsed.small_variant_tsv,
            sv=parsed.sv_result, fusion=parsed.fusion_result,
            viral=parsed.viral_result, isoforms=parsed.isoforms_result)
        parsed.output.parent.mkdir(parents=True, exist_ok=True)
        summary_json = parsed.summary_json or parsed.output.with_suffix(".json")
        summary_json.parent.mkdir(parents=True, exist_ok=True)
        from ..antigen_evidence_pdf import write_antigen_report
        write_antigen_report(parsed.output, inventory)
        summary_json.write_text(json.dumps(inventory, indent=2, sort_keys=True) + "\n")
    except (OSError, ValueError, KeyError, TypeError, json.JSONDecodeError) as error:
        parser_for_program(parser, prog).error(str(error))
    print(parsed.output)
    print(summary_json)
