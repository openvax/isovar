"""Reanalyze source-separated RNA libraries and render complete cohort reports."""

import argparse
import logging

from ..cohort import run_cohort


def run(args=None, prog=None):
    parser = argparse.ArgumentParser(prog=prog or "isovar cohort", description=__doc__)
    parser.add_argument("--manifest", required=True)
    parser.add_argument("--output-dir", required=True)
    parser.add_argument("--resume", action="store_true")
    parser.add_argument("--render-only", action="store_true")
    parser.add_argument("--source-id", action="append")
    parser.add_argument("--no-render", action="store_true", help="Compute only; render all sources together later")
    options = parser.parse_args(args)
    logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s")
    try:
        paths = run_cohort(options.manifest, options.output_dir, resume=options.resume,
                           render_only=options.render_only, source_ids=options.source_id, render=not options.no_render)
    except (OSError, ValueError, KeyError) as error:
        parser.error(str(error))
    for path in paths:
        print(path)
