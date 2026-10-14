"""Prepare a complete event ledger as a checksummed GitHub release attachment."""

import argparse
import os
from pathlib import Path
import re
import shutil
import tempfile

from examples.sid_sv_audit.inventory import digest, read_json, write_json


_MAX_ASSET_BYTES = 2 * 1024 ** 3


def prepare(source, destination, assets, version):
    """Keep full evidence in an asset and produce a separate repository copy."""
    source, destination, assets = map(Path, (source, destination, assets))
    if not re.fullmatch(r"[0-9]+\.[0-9]+\.[0-9]+", version):
        raise ValueError("Expected a stable release version")
    if (destination.resolve().is_relative_to(source.resolve())
            or assets.resolve().is_relative_to(source.resolve())
            or assets.resolve().is_relative_to(destination.resolve())):
        raise ValueError("Report inputs, distribution and assets must use separate locations")
    if destination.exists():
        raise ValueError("Report distribution already exists")
    manifest = source / "files.json"
    manifest_sha = digest(manifest)
    pins = read_json(manifest)
    paths = {str(p.relative_to(source)): p for p in source.rglob("*")
             if p.is_file() and p != source / "files.json"}
    if set(paths) != set(pins):
        raise ValueError("Report files differ from their checksum manifest")
    if "events.json.gz" not in paths:
        raise ValueError("Report is missing its complete event ledger")
    for name, path in paths.items():
        if not path.resolve().is_relative_to(source.resolve()) or digest(path) != pins[name]:
            raise ValueError("Changed or outside report file: " + name)
    ledger = source / "events.json.gz"
    size = ledger.stat().st_size
    if size >= _MAX_ASSET_BYTES:
        raise ValueError("Complete ledger exceeds the GitHub release-asset size limit")
    filename = "sid-neoorf-events-" + version + ".json.gz"
    url = "https://github.com/openvax/isovar/releases/download/v" + version + "/" + filename
    assets.mkdir(parents=True, exist_ok=True)
    asset = assets / filename
    if asset.exists() and digest(asset) != pins["events.json.gz"]:
        raise ValueError("Refusing to replace a different release asset")
    destination.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(prefix=".report-release-", dir=destination.parent) as temporary:
        staged = Path(temporary) / "report"
        staged.mkdir()
        for name, path in paths.items():
            if name == "events.json.gz":
                continue
            target = staged / name
            target.parent.mkdir(parents=True, exist_ok=True)
            if path.suffix == ".md":
                text = re.sub(r"\((?:\.\./)?events\.json\.gz\)", "(" + url + ")", path.read_text())
                if name == "README.md":
                    text += ("\nThe complete ledger is a release attachment. Save it as `events.json.gz`\n"
                             "in this directory to verify every entry in `files.json`. The download\n"
                             "URL, byte count and checksum are recorded in `downloads.json`.\n")
                target.write_text(text)
            else:
                shutil.copyfile(path, target)
        download = dict(filename=filename, url=url, bytes=size, sha256=pins["events.json.gz"])
        write_json(staged / "downloads.json", dict(
            schema="isovar.sid.report_downloads.v1", version=version,
            source_files_sha256=manifest_sha, files={"events.json.gz": download}))
        checksums = {str(p.relative_to(staged)): digest(p) for p in sorted(staged.rglob("*")) if p.is_file()}
        write_json(staged / "files.json", dict(checksums, **{"events.json.gz": pins["events.json.gz"]}))
        # A hard link publishes the fully copied asset without replacing a
        # concurrent writer's file; copying itself uses bounded library buffers.
        with tempfile.NamedTemporaryFile(dir=assets, prefix=".ledger-", delete=False) as handle:
            copied = Path(handle.name)
        try:
            shutil.copyfile(ledger, copied)
            if digest(copied) != pins["events.json.gz"]:
                raise ValueError("Report ledger changed during copying")
            try:
                os.link(copied, asset)
            except FileExistsError:
                if digest(asset) != pins["events.json.gz"]:
                    raise ValueError("Refusing to replace a different release asset")
            # Recheck every source, including the manifest, before publication.
            current_names = {str(p.relative_to(source)) for p in source.rglob("*")
                             if p.is_file() and p != manifest}
            if (digest(manifest) != manifest_sha or current_names != set(paths)
                    or any(digest(p) != pins[n] for n, p in paths.items())):
                raise ValueError("Report changed during preparation")
            if destination.exists():
                raise ValueError("Report distribution already exists")
            staged.rename(destination)
        finally:
            copied.unlink(missing_ok=True)
    return download


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("source", type=Path)
    parser.add_argument("destination", type=Path)
    parser.add_argument("assets", type=Path)
    parser.add_argument("--version", required=True)
    args = parser.parse_args()
    print(prepare(args.source, args.destination, args.assets, args.version))


if __name__ == "__main__":
    main()
