"""Create, verify, back up and replay portable offline audit bundles."""

import argparse
import json
from pathlib import Path

from ..audit_bundle import AuditBundle, AuditBundleError
from ..audit_replay import replay_sv_rna
from ..sid_audit_bundle import bundle_sid_audit


def make_parser(prog="isovar audit-bundle"):
    parser = argparse.ArgumentParser(prog=prog, description=__doc__)
    commands = parser.add_subparsers(dest="command", required=True)
    create = commands.add_parser("create", help="Freeze explicitly selected files from a JSON specification")
    create.add_argument("--spec", required=True, type=Path)
    create.add_argument("--output", required=True, type=Path)
    sid = commands.add_parser("import-sid", help="Freeze existing Sid audit receipts, subsets and replay plans")
    sid.add_argument("--audit", required=True, type=Path)
    sid.add_argument("--output", required=True, type=Path)
    sid.add_argument("--extra-file", action="append", default=[], metavar="NAME=PATH")
    for command in ("inspect", "verify", "pack", "replay"):
        sub = commands.add_parser(command, help={
            "inspect": "List bundle identity, files, replay runs and preserved outcomes",
            "verify": "Check every file's size and checksum",
            "pack": "Create a verified tar.gz backup",
            "replay": "Reconstruct frozen SV views offline and compare full results",
        }[command])
        sub.add_argument("bundle", type=Path)
        if command in ("pack", "replay"):
            sub.add_argument("--output", required=True, type=Path)
        if command == "replay":
            runs = sub.add_mutually_exclusive_group(required=True)
            runs.add_argument("--run-id", action="append")
            runs.add_argument("--all", action="store_true")
            sub.add_argument("--allow-software-change", action="store_true")
            sub.add_argument("--scratch-dir", type=Path)
    unpack = commands.add_parser("unpack", help="Restore and verify a backup in a new directory")
    unpack.add_argument("archive", type=Path)
    unpack.add_argument("--output", required=True, type=Path)
    return parser


def run(args=None, prog=None):
    """Run the installed CLI; emit JSON on stdout and errors with exit status 2."""
    parser = make_parser(prog or "isovar audit-bundle")
    args = parser.parse_args(args)
    try:
        if args.command == "create":
            spec = json.loads(args.spec.read_text())
            if not isinstance(spec, dict) or set(spec) - {"files", "metadata", "expected_sha256"}:
                raise AuditBundleError("Unexpected create specification fields")
            if not isinstance(spec.get("files"), dict) or any(not isinstance(path, str) for path in spec["files"].values()):
                raise AuditBundleError("Specification files must map logical names to local path strings")
            files = {name: args.spec.parent / path for name, path in spec["files"].items()}
            bundle = AuditBundle.create(args.output, files, spec.get("metadata"),
                                        expected_sha256=spec.get("expected_sha256"))
        elif args.command == "import-sid":
            extra = {}
            for value in args.extra_file:
                name, separator, path = value.partition("=")
                if not separator or not name or not path or name in extra:
                    raise AuditBundleError("Use distinct --extra-file NAME=PATH values")
                extra[name] = Path(path)
            bundle = bundle_sid_audit(args.audit, args.output, extra_files=extra)
        elif args.command == "unpack":
            bundle = AuditBundle.unpack(args.archive, args.output)
        else:
            bundle = AuditBundle.open(args.bundle)
        if args.command == "inspect":
            metadata = bundle.manifest["metadata"]
            result = dict(bundle.verify(), software=bundle.manifest["software"],
                          kind=metadata.get("kind"), runs=sorted(metadata.get("sv_rna_runs", {})),
                          outcomes=metadata.get("outcomes", []))
        elif args.command == "pack":
            result = dict(bundle.verify(), archive=str(bundle.pack(args.output)))
        elif args.command == "replay":
            ids = sorted(bundle.manifest["metadata"].get("sv_rna_runs", {})) if args.all else args.run_id
            if not ids:
                raise AuditBundleError("Bundle has no SV RNA replay runs")
            result = []
            for run_id in ids:
                checkpoint = replay_sv_rna(bundle, run_id, args.output,
                    allow_software_change=args.allow_software_change, scratch_dir=args.scratch_dir)
                result.append(dict(run_id=run_id, comparison=checkpoint["comparison"],
                                   result_sha256=checkpoint["result_sha256"]))
        else:
            result = bundle.verify()
        print(json.dumps(result, sort_keys=True))
    except (OSError, ValueError, TypeError, KeyError) as error:
        parser.error(str(error))


if __name__ == "__main__":
    run()
