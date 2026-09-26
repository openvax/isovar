"""Isovar's Sid test reads, from osteosarc's shared OpenVax bundle, openvax-v1.

Isovar's regression reads from the public Sid osteosarcoma data (CC0) are
members of the ``openvax-v1`` bundle named ``isovar/<path>``, after the files
Isovar's tests read (``<path>`` is relative to ``tests/data``). `path` downloads
and verifies the bundle once, into osteosarc's cache, exports Isovar's read
files in their original formats, and returns a local file. Tests need network
access only the first time.

An exported file holds exactly the member's original records, but coordinate
sorted and under the source's full header, so its bytes differ from any older
copy. Reads embedded in JSON fixtures stay in those fixtures; tests check them
against the bundle.

The acquisition helpers below (`open_dataset`, `fetch_metadata`,
`extract_regions`) serve the regional audit builders, which need an osteosarc
metadata snapshot.
"""

from osteosarc.legacy_fixtures import (
    read_json as read_json,
    write_json as write_json,
    sam_digest as sam_digest,
    sam_regions as sam_regions,
    minimal_header as minimal_header,
)

import argparse
from functools import lru_cache
from hashlib import sha256
import json
import os
from pathlib import Path
import shutil
import tempfile


BUNDLE = "openvax-v1"
# The exact published bundle Isovar's tests were checked against (its manifest's SHA-256).
BUNDLE_MANIFEST_SHA256 = "193623c040fa85e9dae5733fe0939358b9b7119da0e6563183bf5a9727f2d7d4"
PREFIX = "isovar/"
# Bump when the export layout changes, so older exports are never reused.
EXPORT_LAYOUT = 1
# Exported file formats by suffix.
_FORMATS = ((".sam.gz", "sam.gz"), (".sam", "sam"), (".bam", "bam"))

# Fetch failures, per cache, so an offline test run fails once per process, not per test.
_failures = {}


class SidDataUnavailable(RuntimeError):
    """openvax-v1 couldn't be fetched, or isn't the bundle Isovar expects."""


def _key(cache):
    return cache if isinstance(cache, (str, Path, type(None))) else id(cache)


def _cache(cache):
    from osteosarc import Cache
    return cache if isinstance(cache, Cache) else Cache(cache)


def bundle(cache=None):
    """The verified ``openvax-v1`` folder, downloaded into osteosarc's cache on first use."""
    from osteosarc import fetch_bundle
    key = _key(cache)
    if key in _failures:
        raise _failures[key]
    try:
        root = Path(fetch_bundle(BUNDLE, cache=cache))
    except Exception as error:
        _failures[key] = SidDataUnavailable(
            "Could not fetch osteosarc's %s test data (%s). Tests need network access once; to run "
            "offline, set OSTEOSARC_CACHE to a cache that already holds it." % (BUNDLE, error))
        raise _failures[key] from error
    if sha256((root / "manifest.json").read_bytes()).hexdigest() != BUNDLE_MANIFEST_SHA256:
        raise SidDataUnavailable(
            "The installed osteosarc publishes a different %s (%s); Isovar expects manifest %s"
            % (BUNDLE, root, BUNDLE_MANIFEST_SHA256))
    return root


@lru_cache(maxsize=4)
def _members(key, cache):
    from osteosarc import list_bundle
    return tuple(sorted(name[len(PREFIX):] for name in list_bundle(bundle(cache)) if name.startswith(PREFIX)))


def members(cache=None):
    """Isovar's member names in the bundle, without the ``isovar/`` prefix.

    A name with ``#`` selects records inside a file that stays in the
    repository: a JSON pointer, or one read.
    """
    return list(_members(_key(cache), cache))


def _format(name):
    return next((fmt for suffix, fmt in _FORMATS if name.endswith(suffix)), None)


def _digests(folder, names):
    return {name: sha256((folder / name).read_bytes()).hexdigest() for name in names}


@lru_cache(maxsize=4)
def _exported(key, cache):
    """Isovar's read files, exported under their original names once per bundle and layout."""
    from osteosarc import export_bundle
    from osteosarc.cache import file_lock
    root = bundle(cache)
    names = [name for name in members(cache) if "#" not in name]
    expected = {name for name in names} | {name + ".bai" for name in names if _format(name) == "bam"}
    exports = _cache(cache).root / "isovar" / ("sid-%s-files-v%d" % (root.name, EXPORT_LAYOUT))
    exports.parent.mkdir(parents=True, exist_ok=True)
    with file_lock(exports.with_name(exports.name + ".lock")):
        if not exports.exists():
            work = Path(tempfile.mkdtemp(dir=exports.parent, prefix=".sid-export-"))
            try:
                for fmt in {_format(name) for name in names}:
                    chosen = [PREFIX + name for name in names if _format(name) == fmt]
                    for member, path in export_bundle(root, work / "raw", members=chosen, format=fmt).items():
                        target = work / "files" / member[len(PREFIX):]
                        target.parent.mkdir(parents=True, exist_ok=True)
                        os.replace(path, target)
                        if fmt == "bam":
                            os.replace(str(path) + ".bai", str(target) + ".bai")
                files = work / "files"
                (files / "digests.json").write_text(json.dumps(_digests(files, sorted(expected)), indent=1) + "\n")
                for path in files.rglob("*"):
                    if path.is_file():
                        path.chmod(0o444)
                os.rename(files, exports)
            finally:
                shutil.rmtree(work, ignore_errors=True)
        _verify_exports(exports, expected)
    return exports


def _verify_exports(exports, expected):
    """Every expected file is present and unchanged since export."""
    try:
        recorded = json.loads((exports / "digests.json").read_text())
        intact = set(recorded) == set(expected) and _digests(exports, sorted(expected)) == recorded
    except (OSError, ValueError):
        intact = False
    if not intact:
        raise SidDataUnavailable("%s doesn't match the exported reads; delete it to export them again" % exports)


def path(name, cache=None):
    """
    Local copy of one of Isovar's Sid read files, exported from openvax-v1.

    Parameters
    ----------
    name : str
        The file's path under ``tests/data``, such as
        ``"osteosarc/bulk_star_t0.sam.gz"``.
    cache : str or osteosarc.Cache, optional
        osteosarc's cache root; by default ``OSTEOSARC_CACHE``, else the shared
        OpenVax cache.

    Returns
    -------
    pathlib.Path
        The file, with its ``.bai`` index beside it for BAM.

    Raises
    ------
    KeyError
        If ``name`` is not one of Isovar's read files. Selections inside files
        (names with ``#``) are exported with `export` instead.
    """
    if "#" in name:
        raise KeyError("%s selects reads inside a file; use sid_data.export" % name)
    if name not in members(cache):
        raise KeyError("%s is not an Isovar read file in %s" % (name, BUNDLE))
    return _exported(_key(cache), cache) / name


def export(name, output, cache=None):
    """
    Write one of Isovar's members, whole files and ``#`` selections alike, to a new file.

    Parameters
    ----------
    name : str
        A member name, as `members` lists it.
    output : str or Path
        The file to write, ending ``.bam`` (written with its ``.bai`` index),
        ``.sam`` or ``.sam.gz``, in an existing folder. Existing files are never
        overwritten.
    cache : str or osteosarc.Cache, optional

    Returns
    -------
    pathlib.Path
    """
    from osteosarc import export_bundle
    output = Path(output)
    fmt = _format(output.name)
    if fmt is None:
        raise ValueError("Export to a .bam, .sam or .sam.gz file, not %s" % output.name)
    if not output.parent.is_dir():
        raise FileNotFoundError(output.parent)
    if name not in members(cache):
        raise KeyError("%s is not an Isovar member of %s" % (name, BUNDLE))
    targets = [output] + ([Path(str(output) + ".bai")] if fmt == "bam" else [])
    for target in targets:
        if target.exists():
            raise FileExistsError(target)
    with tempfile.TemporaryDirectory(dir=output.parent, prefix=".isovar-export-") as work:
        written = export_bundle(bundle(cache), work, members=[PREFIX + name], format=fmt)
        if not written:
            raise KeyError("%s has no reads in %s" % (name, BUNDLE))
        source, = written.values()
        # A hard link never replaces an existing file, even one created meanwhile.
        for suffix, target in zip(("", ".bai"), targets):
            os.link(str(source) + suffix, target)
    return output


def open_dataset(snapshot=None, cache=None, offline=False):
    """Open an explicitly created snapshot; never create one implicitly."""
    name = snapshot or os.environ.get("ISOVAR_SID_SNAPSHOT")
    if not name:
        raise ValueError("Specify --snapshot or ISOVAR_SID_SNAPSHOT (an osteosarc snapshot name)")
    return _open_dataset(name, cache or os.environ.get("ISOVAR_SID_CACHE"), offline)


@lru_cache(maxsize=4)
def _open_dataset(name, cache, offline):
    from osteosarc import Cache, Dataset
    return Dataset.open(name, cache=Cache(cache, offline=offline), offline=offline)


def fetch_metadata(url):
    """Use osteosarc for bounded metadata/index downloads, never full BAMs."""
    from urllib.parse import urlsplit
    if urlsplit(url).path.lower().endswith((".bam", ".cram", ".sam")):
        raise ValueError("Whole alignments are prohibited; use extract_regions")
    dataset = open_dataset()
    for receipt in dataset.manifest["sources"].values():
        if receipt["url"] == url:
            return dataset.cache.path(receipt), receipt
    if "sid-sijbrandij-osteosarc-dataset" in url:
        asset = dataset.file(url)
        path = dataset.download(asset)
        receipt = dict(url=url, sha256=sha256(path.read_bytes()).hexdigest(),
                       size=path.stat().st_size, snapshot_id=dataset.id, asset_id=asset.id)
    else:
        receipt = dataset.cache.fetch(url, max_bytes=256_000_000)
        path = dataset.cache.path(receipt)
    return path, receipt.to_dict() if hasattr(receipt, "to_dict") else receipt


def extract_regions(url, regions, assembly, output=None, *, dataset=None,
                    reference_lengths=None, timeout=600):
    """Delegate header/index/range acquisition to the snapshot's exact Asset."""
    dataset = dataset or open_dataset()
    asset = dataset.file(url)
    subset = dataset.extract_reads(
        asset, sam_regions(regions, assembly, reference_lengths), timeout=timeout)
    if output is not None:
        output = Path(output)
        if output.exists() or Path(str(output) + ".bai").exists():
            raise FileExistsError(output)
        shutil.copyfile(subset.path, output)
        shutil.copyfile(subset.index_path, str(output) + ".bai")
    return subset


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("command", choices=("list", "path", "export"),
                        help="list members; print a read file's local path; export any member")
    parser.add_argument("name", nargs="?", help="A member, as listed; path takes whole read files only")
    parser.add_argument("--output", type=Path, help="export: the .bam, .sam or .sam.gz file to write")
    parser.add_argument("--cache", help="osteosarc cache root")
    args = parser.parse_args(argv)
    if args.command == "list":
        for name in members(args.cache):
            print(name)
    elif args.name is None:
        parser.error("%s requires a member name" % args.command)
    elif args.command == "path":
        if "#" in args.name:
            parser.error("path takes whole read files; use export for %s" % args.name)
        print(path(args.name, args.cache))
    else:
        if args.output is None:
            parser.error("export requires --output")
        export(args.name, args.output, args.cache)


if __name__ == "__main__":
    main()
