"""Reanalyze source-separated RNA libraries and render complete cohort reports."""

import argparse

from ..cohort import run_cohort
from ..logging import configure_cli_logging
from .output_args import add_log_level_arg
from .validation import input_error


def run(args=None, prog=None):
    parser = argparse.ArgumentParser(prog=prog or "isovar cohort", description=__doc__)
    parser.add_argument("--manifest", required=True,
                        help="Cohort manifest JSON: samples, source-separated RNA libraries and tasks "
                             "(see docs/cohort-reports.md)")
    parser.add_argument("--output-dir", required=True,
                        help="Directory for per-source results, receipts and the cohort PDF reports")
    parser.add_argument("--resume", action="store_true", help="Resume analysis only when recorded inputs, settings and code match")
    parser.add_argument("--render-only", action="store_true", help="Render recorded results after checking manifest identity and artifact digests")
    parser.add_argument("--source-id", action="append",
                        help="Compute only this source (repeatable); needs --no-render, since reports "
                             "cover every source")
    parser.add_argument("--no-render", action="store_true", help="Compute only; render all sources together later")
    add_log_level_arg(parser)
    options = parser.parse_args(args)
    if options.render_only and options.resume:
        parser.error("--render-only renders recorded results; it cannot be combined with --resume")
    configure_cli_logging(options.log_level)
    try:
        paths = run_cohort(options.manifest, options.output_dir, resume=options.resume,
                           render_only=options.render_only, source_ids=options.source_id, render=not options.no_render)
    except KeyError as error:
        input_error(parser, "Missing manifest entry %s" % error)
    except (OSError, ValueError) as error:
        input_error(parser, error)
    for path in paths:
        print(path)
