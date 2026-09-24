"""Explicit acquisition of the pinned osteosarc original-read regression subset.

Run ``python -m isovar.osteosarc_data --help``. Importing this module does no
I/O; ordinary Isovar tests use their checked-in fixtures without acquisition.
"""

import argparse
from hashlib import sha256
import json
from pathlib import Path
import re
import shutil
import subprocess
from urllib.parse import urlsplit


DEFAULT_MANIFEST = Path(__file__).with_name("data") / "osteosarc-test-data-v1.json"


def load_manifest(path=DEFAULT_MANIFEST):
    """Read the consumer-owned dataset identity and immutable asset pins."""
    manifest = json.loads(Path(path).read_text())
    if (manifest.get("schema_version") != 1 or manifest.get("dataset") != "osteosarc"
            or not manifest.get("data_version") or not manifest.get("assets")):
        raise ValueError("Unsupported osteosarc fixture manifest")
    names = set()
    for asset in manifest["assets"]:
        name = asset["filename"]
        if (not isinstance(name, str) or name in ("", ".", "..", "manifest.json")
                or any(c in name for c in "/\\:") or name in names):
            raise ValueError("Asset filenames must be unique plain basenames")
        names.add(name)
        if (not re.fullmatch(r"[0-9a-f]{64}", asset["sha256"])
                or type(asset["size_bytes"]) is not int or asset["size_bytes"] < 0
                or urlsplit(asset["url"]).scheme != "https"
                or Path(urlsplit(asset["url"]).path).suffixes != Path(name).suffixes):
            raise ValueError("Asset needs a SHA-256, byte size and HTTPS URL with matching suffixes")
    cases = manifest.get("cases", [])
    if (not cases or len({c["case_id"] for c in cases}) != len(cases)
            or any(c["bam"] not in names or c["bam"] + ".bai" not in names for c in cases)):
        raise ValueError("Cases must be unique and include their BAM and index assets")
    return manifest


def _verify_file(path, asset):
    if path.is_symlink() or not path.is_file():
        raise ValueError("Missing or symlinked fixture object: %s" % path)
    if path.stat().st_size != asset["size_bytes"] or sha256(path.read_bytes()).hexdigest() != asset["sha256"]:
        raise ValueError("Fixture checksum mismatch: %s" % path)


def verify_dataset(output, manifest_path=DEFAULT_MANIFEST):
    """Verify a complete export without network access, repair or cache writes."""
    manifest, output = load_manifest(manifest_path), Path(output)
    if output.is_symlink() or not output.is_dir():
        raise ValueError("Dataset must be a real directory")
    if {p.name for p in output.iterdir()} != {a["filename"] for a in manifest["assets"]} | {"manifest.json"}:
        raise ValueError("Dataset has missing or unexpected files")
    pin = output / "manifest.json"
    if pin.is_symlink() or json.loads(pin.read_text()) != manifest:
        raise ValueError("Export manifest does not match the requested dataset")
    for asset in manifest["assets"]:
        _verify_file(output / asset["filename"], asset)
    return output


def acquire_dataset(output=None, *, manifest_path=DEFAULT_MANIFEST, cache_root=None,
                    import_corpus=None, offline=False, repair_cache=False):
    """Import or fetch pinned objects, optionally exporting a new offline subset.

    Valid shared objects are adopted with osteosarc import receipts without
    rewriting their bytes. Invalid objects fail unless online repair is
    explicitly requested. Existing exports are only verified, never replaced.
    ``import_corpus`` supplies original filenames from an existing checkout;
    it does not change the manifest's selection or allele definitions.
    """
    if offline and repair_cache:
        raise ValueError("Offline mode cannot repair the cache")
    manifest = load_manifest(manifest_path)
    output = Path(output) if output is not None else None
    if output is not None and (output.exists() or output.is_symlink()):
        return verify_dataset(output, manifest_path)
    from osteosarc import Cache, OsteosarcError
    cache = Cache(cache_root, offline=offline)
    paths = {}
    for asset in manifest["assets"]:
        cached = cache.objects / (asset["sha256"] + "".join(Path(asset["filename"]).suffixes))
        source = cached if cached.exists() or cached.is_symlink() else (
            Path(import_corpus) / asset["filename"] if import_corpus is not None else None)
        repair = False
        if source is not None:
            try:
                _verify_file(source, asset)
            except ValueError:
                if source != cached or not repair_cache or cached.is_symlink():
                    raise
                repair = True
        expected = dict(sha256=asset["sha256"], size=asset["size_bytes"])
        try:
            receipt = (cache.import_file(source, asset["url"], **expected) if source is not None and not repair
                       else cache.fetch(asset["url"], refresh=repair_cache, **expected))
            paths[asset["filename"]] = cache.path(receipt)
        except OsteosarcError as error:
            raise ValueError(str(error)) from error
    if output is None:
        return paths
    # Resolve every input first; a missing/corrupt cache cannot publish an export.
    # Exclusive creation protects existing directories, including concurrent exports.
    output.mkdir(parents=True)
    for name, path in paths.items():
        with path.open("rb") as source, (output / name).open("xb") as target:
            shutil.copyfileobj(source, target)
    with (output / "manifest.json").open("x") as handle:
        json.dump(manifest, handle, indent=2, sort_keys=True)
        handle.write("\n")
    return verify_dataset(output, manifest_path)


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--manifest", type=Path, default=DEFAULT_MANIFEST)
    parser.add_argument("--cache-root", type=Path, help="Shared cache root; defaults to osteosarc/OpenVax settings")
    parser.add_argument("--import-corpus", type=Path, help="Import pinned files from an existing fixture directory")
    parser.add_argument("--output", type=Path, help="Export a new directory, or verify an existing identical export")
    parser.add_argument("--offline", action="store_true", help="Use local objects only")
    parser.add_argument("--repair-cache", action="store_true", help="Explicitly refetch corrupt cache objects")
    parser.add_argument("--verify-only", action="store_true", help="Verify --output without touching the cache")
    args = parser.parse_args(argv)
    if args.verify_only and (args.output is None or args.import_corpus or args.repair_cache):
        parser.error("--verify-only requires --output and cannot import or repair")
    try:
        manifest = load_manifest(args.manifest)
        result = (verify_dataset(args.output, args.manifest) if args.verify_only else acquire_dataset(
            args.output, manifest_path=args.manifest, cache_root=args.cache_root, import_corpus=args.import_corpus,
            offline=args.offline, repair_cache=args.repair_cache))
    except (OSError, ValueError, RuntimeError, subprocess.SubprocessError) as error:
        parser.exit(1, "Osteosarc fixtures unavailable: %s\n" % error)
    print(json.dumps(dict(dataset=manifest["dataset"], data_version=manifest["data_version"],
                          cases=len(manifest["cases"]), assets=len(manifest["assets"]),
                          output=str(result) if isinstance(result, Path) else None), indent=2))


if __name__ == "__main__":
    main()
