"""Run a reproducible, source-separated nominated-SV RNA audit."""

from ..logging import configure_cli_logging, get_logger
from .help import HelpArgumentParser
from .input_args import add_library_cache_arg, add_manifest_arg
from .output_args import add_log_level_arg, add_output_dir_arg
from .validation import input_error

logger = get_logger(__name__)


def make_parser(prog="isovar sv-audit"):
    parser = HelpArgumentParser(prog=prog, description=__doc__)
    add_manifest_arg(
        parser,
        required=True,
        help="Source-pinned audit manifest JSON: files with checksums, sources, events and "
             "views (see docs/sv-rna-audit.md)",
    )
    add_output_dir_arg(
        parser,
        required=True,
        help="Directory for per-run snapshots, receipts and report.json/report.tsv",
    )
    parser.add_argument("--resume", action="store_true",
                        help="Reuse snapshots whose implementation, inputs and parameters match")
    parser.add_argument("--event-id", action="append", help="Select event IDs; default: every manifest event")
    add_library_cache_arg(
        parser,
        dest="library_cache_dir",
        aliases=("--library-cache-dir",),
        help="Retain verified original-RNA indexes for reuse across audit invocations",
    )
    add_log_level_arg(parser)
    return parser


def run(args=None, prog=None):
    parser = make_parser(prog or "isovar sv-audit")
    options = parser.parse_args(args)
    configure_cli_logging(options.log_level)
    try:
        from ..sv_rna_audit import run_sv_rna_audit

        report = run_sv_rna_audit(options.manifest, options.output_dir, resume=options.resume,
                                  event_ids=options.event_id, library_cache_dir=options.library_cache_dir)
        logger.info("Assessed %d event/source/view runs; report in %s", len(report["runs"]), options.output_dir)
    except KeyError as error:
        input_error(parser, "Missing manifest entry %s" % error)
    except (OSError, ValueError) as error:
        input_error(parser, error)
