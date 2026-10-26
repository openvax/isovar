"""Run a reproducible, source-separated nominated-SV RNA audit."""

import argparse

from ..sv_rna_audit import run_sv_rna_audit


def run(args=None, prog=None):
    parser = argparse.ArgumentParser(prog=prog or "isovar sv-audit", description=__doc__)
    parser.add_argument("--manifest", required=True)
    parser.add_argument("--output-dir", required=True)
    parser.add_argument("--resume", action="store_true")
    parser.add_argument("--event-id", action="append", help="Select event IDs; default: every manifest event")
    options = parser.parse_args(args)
    try:
        run_sv_rna_audit(options.manifest, options.output_dir,
                                resume=options.resume, event_ids=options.event_id)
    except (OSError, ValueError, KeyError) as error:
        parser.error(str(error))
