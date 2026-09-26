"""Isovar's Sid test reads, from osteosarc's shared OpenVax bundle, openvax-v1.

Isovar's regression reads from the public Sid osteosarcoma data (CC0) are
members of the ``openvax-v1`` bundle named ``isovar/<path>``, after the files
Isovar's tests read (``<path>`` is relative to ``tests/data``). `path` returns
one as a local file, through osteosarc's ``bundle_file``: the bundle is
downloaded and verified once, and each member is exported once into
osteosarc's cache, read-only, and reused offline. Tests need network access
only the first time.

An exported file holds exactly the member's original records, but coordinate
sorted and under the source's full header, and is named ``<member>.<format>``.
Reads embedded in JSON fixtures stay in those fixtures; tests check them
against the bundle.

The acquisition helpers below (`open_dataset`, `fetch_metadata`,
`extract_regions`) serve the regional audit builders, which need an osteosarc
metadata snapshot.
"""

import argparse
from functools import lru_cache
import gzip
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
# Exported file formats by suffix.
_FORMATS = ((".sam.gz", "sam.gz"), (".sam", "sam"), (".bam", "bam"))

# Fetch failures, per cache, so an offline test run fails once per process, not per test.
_failures = {}


class SidDataUnavailable(RuntimeError):
    """openvax-v1 couldn't be fetched, or isn't the bundle Isovar expects."""


def _key(cache):
    return cache if isinstance(cache, (str, Path, type(None))) else id(cache)


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


def _file(name, fmt, cache):
    from osteosarc import bundle_file
    return Path(bundle_file(bundle(cache), PREFIX + name, format=fmt, cache=cache))


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
        A read-only file in the original format, named ``<member>.<format>``,
        with its ``.bai`` index beside it for BAM.

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
    return _file(name, _format(name), cache)


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
    source = _file(name, fmt, cache)
    for suffix, target in zip(("", ".bai"), targets):
        # "x" never replaces a file, even one created meanwhile.
        with open(str(source) + suffix, "rb") as reader, open(target, "xb") as writer:
            shutil.copyfileobj(reader, writer)
    return output


def read_json(path):
    """JSON from a file, gzip-compressed when its name ends ``.gz``."""
    with gzip.open(path, "rt") if str(path).endswith(".gz") else open(path) as handle:
        return json.load(handle)


def write_json(path, value):
    """Publish deterministic gzip JSON atomically, never overwriting a file."""
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(dir=path.parent, prefix=".sid-json-") as temporary:
        staged = Path(temporary) / "data.gz"
        with staged.open("wb") as raw:
            with gzip.GzipFile(filename="", mode="wb", fileobj=raw, mtime=0) as handle:
                handle.write(json.dumps(value, sort_keys=True, separators=(",", ":")).encode())
        if staged.stat().st_size > 64 * 1024 * 1024:
            raise ValueError("Fixture exceeds its 64 MiB size budget")
        os.link(staged, path)


def sam_digest(line):
    """SHA-256 of one SAM record's text."""
    return sha256(line.encode("ascii")).hexdigest()


def sam_regions(regions, assembly, reference_lengths=None):
    """Convert explicit 1-based inclusive SAM intervals (``contig:start-end``) to osteosarc Regions."""
    from osteosarc import Region
    result = []
    for region in regions:
        contig, span = region.rsplit(":", 1)
        start, end = map(int, span.split("-"))
        result.append(Region(contig, start - 1, end, assembly,
                             reference_length=(reference_lengths or {}).get(contig)))
    if not result:
        raise ValueError("An explicit nonempty list of locus intervals is required")
    return result


def minimal_header(header, records):
    """
    The part of a SAM header that decodes ``records``: their references
    (including SA targets), read groups, and program chains.
    """
    references, groups, programs = set(), set(), set()
    for line in records:
        fields = line.split("\t")
        references.update(v for v in (fields[2], fields[6]) if v not in ("*", "="))
        groups.update(v[5:] for v in fields[11:] if v.startswith("RG:Z:"))
        programs.update(v[5:] for v in fields[11:] if v.startswith("PG:Z:"))
        for tag in fields[11:]:
            if tag.startswith("SA:Z:"):
                references.update(row.split(",", 1)[0] for row in tag[5:].split(";") if row)
    sq = [row for row in header.get("SQ", []) if row["SN"] in references]
    if {row["SN"] for row in sq} != references:
        raise ValueError("Source header is missing a required reference sequence")
    # Some source records have RG tags absent from their header. Keep that
    # condition rather than inventing read-group metadata.
    result = {"HD": {"VN": "1.6", "SO": "unknown"}, "SQ": sq}
    rg = [row for row in header.get("RG", []) if row["ID"] in groups]
    if rg:
        result["RG"] = rg
    programs.update(row["PG"] for row in rg if "PG" in row)
    by_id = {row["ID"]: row for row in header.get("PG", [])}
    pending = list(programs)
    while pending:
        previous = by_id.get(pending.pop(), {}).get("PP")
        if previous and previous not in programs:
            programs.add(previous)
            pending.append(previous)
    pg = [row for row in header.get("PG", []) if row["ID"] in programs]
    if pg:
        result["PG"] = pg
    return result


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
