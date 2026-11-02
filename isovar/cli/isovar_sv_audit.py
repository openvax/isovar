"""Run a reproducible, source-separated nominated-SV RNA audit."""

import argparse

from .output_args import add_log_level_arg
from ..logging import configure_cli_logging, get_logger
from ..sv_rna_audit import run_sv_rna_audit
from .validation import input_error

logger = get_logger(__name__)


def run(args=None, prog=None):
    parser = argparse.ArgumentParser(prog=prog or "isovar sv-audit", description=__doc__)
    parser.add_argument("--manifest", required=True,
                        help="Source-pinned audit manifest JSON: files with checksums, sources, events and "
                             "views (see docs/sv-rna-audit.md)")
    parser.add_argument("--output-dir", required=True,
                        help="Directory for per-run snapshots, receipts and report.json/report.tsv")
    parser.add_argument("--resume", action="store_true",
                        help="Reuse snapshots whose implementation, inputs and parameters match")
    parser.add_argument("--event-id", action="append", help="Select event IDs; default: every manifest event")
    parser.add_argument("--library-cache-dir", help="Retain verified original-RNA indexes for reuse across audit invocations")
    add_log_level_arg(parser)
    options = parser.parse_args(args)
    configure_cli_logging(options.log_level)
    try:
        report = run_sv_rna_audit(options.manifest, options.output_dir, resume=options.resume,
                                  event_ids=options.event_id, library_cache_dir=options.library_cache_dir)
        logger.info("Assessed %d event/source/view runs; report in %s", len(report["runs"]), options.output_dir)
    except KeyError as error:
        input_error(parser, "Missing manifest entry %s" % error)
    except (OSError, ValueError) as error:
        input_error(parser, error)
