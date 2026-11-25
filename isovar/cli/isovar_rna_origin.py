"""Assess nuclear/mitochondrial RNA origin relative to explicit competing paths."""

import json
from pathlib import Path

from ..default_parameters import RNA_ORIGIN_DEFAULTS
from .help import HelpArgumentParser
from .input_args import add_genome_fasta_arg, add_sample_id_arg, add_source_arg
from .rna_args import add_rna_path_args, rna_path_context
from .validation import input_error
from .output_args import add_force_arg, add_log_level_arg, add_output_file_arg
from ..logging import configure_cli_logging


def make_parser(prog="isovar rna-origin"):
    parser = HelpArgumentParser(prog=prog, description=__doc__)
    add_rna_path_args(parser)
    add_genome_fasta_arg(parser, help="Reference FASTA when needed to read CRAM")
    parser.add_argument("--hypotheses", required=True, help="Origin-panel JSON; see docs/rna-origin.md")
    add_sample_id_arg(parser, required=True, help="Sample providing the RNA")
    add_source_arg(
        parser,
        required=True,
        dest="source_id",
        aliases=("--source-id",),
        help="Stable original RNA library identity",
    )
    add_output_file_arg(parser, required=True, help="New JSON output file with heuristic origin evidence")
    add_force_arg(parser, help="Overwrite an existing --output file")
    for name in ("min_base_quality", "min_context_bases", "min_edit_margin", "max_edit_fraction"):
        parser.add_argument("--" + name.replace("_", "-"),
                            type=float if name == "max_edit_fraction" else int,
                            default=RNA_ORIGIN_DEFAULTS[name], help="Origin heuristic threshold (default %(default)s)")
    add_log_level_arg(parser)
    return parser


def run(args=None, prog=None):
    parser = make_parser(prog or "isovar rna-origin")
    options = parser.parse_args(args)
    configure_cli_logging(options.log_level)
    try:
        from ..json_files import publish_json
        from ..rna_origin import RNAOriginHypothesis, RNAOriginPolicy, assess_rna_origin

        output = Path(options.output)
        if output.exists() and not options.force:
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
