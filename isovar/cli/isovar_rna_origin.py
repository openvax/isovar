"""Assess nuclear/mitochondrial RNA origin relative to explicit competing paths."""

import argparse
import json
from pathlib import Path

from .rna_args import add_rna_path_args, rna_path_context
from .validation import input_error
from .output_args import add_log_level_arg
from ..logging import configure_cli_logging
from ..audit_bundle import publish_json
from ..rna_origin import RNAOriginHypothesis, RNAOriginPolicy, assess_rna_origin


def make_parser(prog="isovar rna-origin"):
    parser = argparse.ArgumentParser(prog=prog, description=__doc__)
    add_rna_path_args(parser)
    parser.add_argument("--genome-fasta", help="Reference FASTA when needed to read CRAM")
    parser.add_argument("--hypotheses", required=True, help="Origin-panel JSON; see docs/rna-origin.md")
    parser.add_argument("--sample-id", required=True, help="Sample providing the RNA")
    parser.add_argument("--source-id", required=True, help="Stable original RNA library identity")
    parser.add_argument("--output", required=True, help="New JSON file with heuristic origin evidence")
    defaults = RNAOriginPolicy()
    for name in ("min_base_quality", "min_context_bases", "min_edit_margin", "max_edit_fraction"):
        parser.add_argument("--" + name.replace("_", "-"),
                            type=float if name == "max_edit_fraction" else int,
                            default=getattr(defaults, name), help="Origin heuristic threshold (default %(default)s)")
    add_log_level_arg(parser)
    return parser


def run(args=None, prog=None):
    parser = make_parser(prog or "isovar rna-origin")
    options = parser.parse_args(args)
    configure_cli_logging(options.log_level)
    try:
        output = Path(options.output)
        if output.exists():
            raise ValueError("Output already exists: %s" % output)
        request = json.loads(Path(options.hypotheses).read_text())
        if (not isinstance(request, dict) or set(request) - {"schema", "hypotheses", "reference_name"}
                or request.get("schema") != "isovar.rna_origin_panel.v1"):
            raise ValueError("Expected an isovar.rna_origin_panel.v1 panel")
        hypotheses = [RNAOriginHypothesis(**row) for row in request["hypotheses"]]
        policy = RNAOriginPolicy(**{name: getattr(options, name) for name in
                                   ("min_base_quality", "min_context_bases", "min_edit_margin", "max_edit_fraction")})
        reference_name = request.get("reference_name")
        with rna_path_context(options, reference_name=reference_name) as (bam, collector, library):
            result = assess_rna_origin(bam, hypotheses, sample_id=options.sample_id, source=options.source_id,
                                      policy=policy, read_collector=collector, library_index=library,
                                      reference_name=reference_name)
        publish_json(output, result)
    except (OSError, ValueError, TypeError, KeyError) as error:
        input_error(parser, error)
