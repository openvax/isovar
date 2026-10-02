"""Run with python -m examples.sid_sv_audit --help."""

import argparse
from pathlib import Path

from . import acquisition, inventory, reconstruct, references, report


def main():
    parser = argparse.ArgumentParser(description="Catalogue-wide Sid SV RNA reconstruction audit")
    parser.add_argument("command", choices=("inventory", "references", "headers", "acquire", "run", "report"))
    parser.add_argument("directory", type=Path)
    parser.add_argument("--snapshot")
    parser.add_argument("--cache", required=True)
    parser.add_argument("--metadata-cache")
    parser.add_argument("--cohort", choices=("tumor_candidate", "blood_control", "unresolved_specimen"),
                        default="tumor_candidate")
    parser.add_argument("--workers", type=int, default=2)
    parser.add_argument("--source-id")
    parser.add_argument("--max-seconds", type=int, default=180)
    parser.add_argument("--output", type=Path)
    parser.add_argument("--require-complete", action="store_true")
    args = parser.parse_args()
    if not 1 <= args.workers <= 8:
        parser.error("Use 1 through 8 workers")
    if args.command == "inventory":
        if not args.snapshot:
            parser.error("inventory needs an explicit --snapshot")
        result = inventory.freeze(args.directory, args.snapshot, args.cache, args.metadata_cache or args.cache)
        print(inventory.summary(result))
    elif args.command == "references":
        result = references.build(args.directory)
        print("Pinned %d transcript models" % len(result["models"]))
    elif args.command == "headers":
        acquisition.survey(args.directory, args.cache, args.cohort, args.workers)
    elif args.command == "acquire":
        acquisition.acquire(args.directory, args.cache, args.cohort, args.workers, args.source_id)
    elif args.command == "run":
        reconstruct.run(args.directory, args.cohort, args.workers, args.source_id, args.max_seconds)
    elif args.command == "report":
        if not args.output:
            parser.error("report needs --output (a new directory)")
        print(report.report(args.directory, args.output, args.cohort, args.require_complete))


if __name__ == "__main__":
    main()
