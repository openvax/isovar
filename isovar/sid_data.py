"""Isovar's Sid test reads, packaged with Isovar and made with osteosarc.

Isovar's regression reads from the public Sid osteosarcoma data (CC0) ship
with Isovar in ``isovar/data/sid-reads``, so its tests never need network
access. That folder is an osteosarc fixture bundle of Isovar's members of
osteosarc's shared OpenVax bundle, ``openvax-v1``. Legacy members are named
``isovar/<path>``, after the files Isovar's tests read (``<path>`` is relative
to ``tests/data``). `path` returns one as a local file, through osteosarc's
``bundle_file``, which exports each member once into osteosarc's cache,
read-only, and reuses it.

The bundle also holds the shared T2 small-variant panel under its original
member names. These are exported with ``osteosarc.bundle_file(bundle(), name)``.

An exported file holds exactly the member's original records, but coordinate
sorted and under the source's full header, and is named ``<member>.<format>``.
JSON fixtures refer to these records by member name and SAM digest. `read_fixture`
restores them without changing record order, multiplicity or adjacent metadata.

`build` makes the packaged bundle again with osteosarc, from a shared bundle's
records or from Sid's original BAMs, and `check` compares such a build with
the packaged bundle.

The acquisition helpers below (`open_dataset`, `fetch_metadata`,
`extract_regions`) serve the regional audit builders, which need an osteosarc
metadata snapshot.
"""

import argparse
from collections import Counter
from functools import lru_cache
import gzip
from hashlib import sha256
import json
import os
from pathlib import Path
import shutil
import sys
import tempfile


BUNDLE = "openvax-v1"
# The shared bundle Isovar's members come from (its manifest's SHA-256).
BUNDLE_MANIFEST_SHA256 = "193623c040fa85e9dae5733fe0939358b9b7119da0e6563183bf5a9727f2d7d4"
# Isovar's own bundle of those members, packaged with Isovar, and its manifest's SHA-256.
PACKAGED = Path(__file__).parent / "data" / "sid-reads"
PACKAGED_MANIFEST_SHA256 = "0e299b970dfaf9df8ca42300db482dd9058642c307ccfca9cca6bc08b078a8c5"
PREFIX = "isovar/"
# Published T2 short-read selection, retained under its upstream member names.
PANEL_SOURCE = "25.03.23.rna.ucla.2025.01.resection.tcga.d32.protocolAligned.sorted"
# Exported file formats by suffix.
_FORMATS = ((".sam.gz", "sam.gz"), (".sam", "sam"), (".bam", "bam"))


class SidDataUnavailable(RuntimeError):
    """Isovar's packaged test reads are missing or changed, or can't be made again."""


def _digest(path):
    return sha256(Path(path).read_bytes()).hexdigest()


def bundle():
    """The packaged bundle's folder, verified once against its pinned manifest hash."""
    _verified(PACKAGED)
    return PACKAGED


@lru_cache(maxsize=2)
def _verified(folder):
    from osteosarc import verify_bundle
    try:
        return verify_bundle(folder, sha256=PACKAGED_MANIFEST_SHA256)
    except (OSError, ValueError) as error:  # osteosarc's IntegrityError is a ValueError.
        raise SidDataUnavailable(
            "Isovar's packaged Sid test reads (%s) are missing or damaged (%s): reinstall Isovar, "
            "or make them again with `python -m isovar.sid_data build`" % (folder, error)) from error


def members():
    """Isovar's member names in the bundle, without the ``isovar/`` prefix.

    A name with ``#`` selects records inside a file that stays in the
    repository: a JSON pointer, or one read.
    """
    return sorted(name[len(PREFIX):] for name in _verified(bundle())["members"] if name.startswith(PREFIX))


def _format(name):
    return next((fmt for suffix, fmt in _FORMATS if name.endswith(suffix)), None)


def _file(name, fmt, cache):
    from osteosarc import bundle_file
    return Path(bundle_file(bundle(), PREFIX + name, format=fmt, cache=cache))


def path(name, cache=None):
    """
    Local copy of one of Isovar's Sid read files, exported from the packaged bundle.

    Parameters
    ----------
    name : str
        The file's path under ``tests/data``, such as
        ``"osteosarc/bulk_star_t0.sam.gz"``.
    cache : str or osteosarc.Cache, optional
        osteosarc's cache root, where exports are kept; by default
        ``OSTEOSARC_CACHE``, else the shared OpenVax cache.

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
    if name not in members():
        raise KeyError("%s is not one of Isovar's Sid read files" % name)
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
    if name not in members():
        raise KeyError("%s is not one of Isovar's Sid test read members" % name)
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


def recipe(shared):
    """
    Isovar's fixtures and the T2 small-variant panel, without reselecting reads.

    Keep all small-variant target definitions, including reference-specific
    targets without a T2 member and unresolved catalog alleles. Their absence
    must remain distinguishable from an observed lack of alternate reads.

    Parameters
    ----------
    shared : dict
        A bundle's recipe, such as openvax-v1's ``recipe.json``.

    Returns
    -------
    dict
        A recipe named ``isovar``, which osteosarc makes Isovar's bundle from.
    """
    panel_targets = {name for name, target in shared["targets"].items()
                     if target["kind"] == "small_variant"
                     or (target["kind"] == "unresolved" and target.get("label") == "current")}
    members = {name: member for name, member in shared["members"].items()
               if name.startswith(PREFIX)
               or (member["source"] == PANEL_SOURCE and member["target"] in panel_targets)}
    sources = {member["source"] for member in members.values()}
    targets = {member["target"] for member in members.values()} | panel_targets
    result = dict(shared, id="isovar", members=members,
                  sources={name: s for name, s in shared["sources"].items() if name in sources},
                  targets={name: t for name, t in shared["targets"].items() if name in targets})
    if "aliases" in shared:
        # Old target names, by the target each now names.
        result["aliases"] = {old: new for old, new in shared["aliases"].items() if new in targets}
    return result


def build(destination, shared=BUNDLE, *, sid=False, cache=None):
    """
    Make Isovar's bundle of Sid test reads again with osteosarc, in a new folder.

    Parameters
    ----------
    destination : str or Path
        A new folder.
    shared : str or Path
        The bundle to take Isovar's members and their recipe from: a published
        bundle's name or a bundle folder (such as `PACKAGED`), as osteosarc
        reads them. openvax-v1 must be the one Isovar pins,
        `BUNDLE_MANIFEST_SHA256`.
    sid : bool
        Acquire every record again from Sid's original BAMs, reading their
        indexed windows over the network, instead of taking the records from
        ``shared``. This needs the osteosarc metadata snapshot the recipe was
        made from.
    cache : str or osteosarc.Cache, optional
        osteosarc's cache root.

    Returns
    -------
    dict
        The new bundle's manifest. To package the bundle, put it in place of
        `PACKAGED` and pin its ``manifest.json``'s SHA-256 as
        `PACKAGED_MANIFEST_SHA256`.
    """
    from osteosarc import generate_bundle, verify_bundle
    from osteosarc.shared import bundle_folder
    folder = Path(bundle_folder(shared, cache=cache))
    whole = read_json(folder / "recipe.json")
    # openvax-v1 must be the one Isovar pins, whatever folder it's in.
    manifest = verify_bundle(folder, sha256=BUNDLE_MANIFEST_SHA256 if whole["id"] == BUNDLE else None)
    mine = recipe(whole)
    if sid:
        from osteosarc import Dataset
        snapshot = whole["snapshot"]
        try:
            # Each source's index comes from the metadata snapshot the recipe was made from.
            dataset, sources = Dataset.open(snapshot["id"], cache=cache, offline=False), None
        except FileNotFoundError as error:
            raise SidDataUnavailable(
                "Acquiring the reads again needs osteosarc's metadata snapshot %s (%s), which %s was made "
                "from (%s). A new snapshot won't do: the recipe pins each BAM to that one."
                % (snapshot["name"], snapshot["id"], whole["id"], error)) from error
    else:
        # Sources whose members are all omitted or unresolved hold no records, and aren't read.
        dataset, sources = None, {name: folder / manifest["sources"][name]["bam"]
                                  for name in mine["sources"] if name in manifest["sources"]}
    return generate_bundle(mine, destination, sources=sources, cache=cache, dataset=dataset,
                           parent=dict(bundle=whole["id"], manifest_sha256=_digest(folder / "manifest.json")))


def check(shared=BUNDLE, *, sid=False, cache=None):
    """
    Make Isovar's bundle again with `build`, and compare it with the packaged one.

    Parameters
    ----------
    shared, sid, cache
        As for `build`.

    Returns
    -------
    list of str
        How the new bundle's recipe, members or sources differ from the
        packaged bundle's: empty when osteosarc makes the same records, under
        the same exported headers, for the same members. How records were
        acquired may differ: the members' ``acquisition_status``, the
        acquisition receipts, and the archived original headers (a bundle
        carved from another keeps its exported headers as the originals).
    """
    packaged = _verified(bundle())
    with tempfile.TemporaryDirectory(prefix="isovar-sid-reads-") as work:
        built = build(Path(work) / "bundle", shared, sid=sid, cache=cache)

    def held(manifest, part, name):
        entry = manifest[part].get(name)
        return entry if part == "sources" or entry is None else {
            key: value for key, value in entry.items() if key != "acquisition_status"}
    problems = ["recipe differs"] if built["recipe_sha256"] != packaged["recipe_sha256"] else []
    for part in ("members", "sources"):
        for name in sorted(built[part].keys() | packaged[part].keys()):
            if held(built, part, name) != held(packaged, part, name):
                problems.append("%s %s differs" % (part[:-1], name))
    return problems


def read_json(path):
    """JSON from a file, gzip-compressed when its name ends ``.gz``."""
    with gzip.open(path, "rt") if str(path).endswith(".gz") else open(path) as handle:
        return json.load(handle)


def records(name, cache=None):
    """Return one member's original SAM text, keyed by its SHA-256 digest.

    Accepts whole-file and ``#`` member names. The packaged bundle and the
    exported records are verified before use; no network access is needed.
    Duplicate records share a digest, while the member's multiplicities are
    checked before constructing the lookup.
    """
    from osteosarc.records import read_records

    manifest = _verified(bundle())
    member = manifest["members"][PREFIX + name]
    rows = list(read_records(_file(name, "bam", cache)))
    # The bundle uses lossless BAM identities, which are deliberately distinct
    # from the SAM-text digests used by the JSON fixtures.
    counts = Counter(row.digest for row in rows)
    if counts != member["records"]:
        raise SidDataUnavailable("Exported Sid records differ from the packaged member: %s" % name)
    lines = [row.read.to_string() for row in rows]
    return {sam_digest(line): line for line in lines}


def restore_records(value, cache=None):
    """Decode explicit Sid record references in a JSON value, leaving it unchanged.

    A record reference is exactly ``{"sid_member": name, "sam_sha256": digest}``.
    Only that member may supply the record. Other dictionaries, lists, scalar
    values and their order are preserved in the returned copy. A missing member
    or digest raises ``KeyError``; malformed references raise ``ValueError``.
    """
    by_member = {}

    def restore(item):
        if isinstance(item, list):
            return [restore(child) for child in item]
        if isinstance(item, dict):
            if "sid_member" in item:
                if (set(item) != {"sid_member", "sam_sha256"}
                        or not isinstance(item["sid_member"], str)
                        or not isinstance(item["sam_sha256"], str)):
                    raise ValueError("Invalid Sid record reference: %r" % item)
                name = item["sid_member"]
                if name not in by_member:
                    by_member[name] = records(name, cache=cache)
                return by_member[name][item["sam_sha256"]]
            return {key: restore(child) for key, child in item.items()}
        return item

    return restore(value)


def read_fixture(path, cache=None):
    """Read plain or gzip JSON and restore its explicit Sid record references."""
    return restore_records(read_json(path), cache=cache)


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
    parser.add_argument("command", choices=("list", "path", "export", "build", "check"),
                        help="list members; print a read file's local path; export any member; make the "
                             "packaged bundle again in a new folder; compare a new build with it")
    parser.add_argument("name", nargs="?", help="A member, as listed; path takes whole read files only")
    parser.add_argument("--output", type=Path,
                        help="export: the .bam, .sam or .sam.gz file to write; build: a new folder")
    parser.add_argument("--from", dest="shared", default=BUNDLE,
                        help="build, check: the bundle to take Isovar's members from, a published "
                             "bundle's name or a folder's path (default: %(default)s)")
    parser.add_argument("--sid", action="store_true",
                        help="build, check: acquire every record again from Sid's original BAMs (slow; "
                             "needs network access and the osteosarc snapshot the recipe names)")
    parser.add_argument("--cache", help="osteosarc cache root")
    args = parser.parse_args(argv)
    if args.command == "list":
        for name in members():
            print(name)
    elif args.command == "build":
        if args.output is None:
            parser.error("build requires --output")
        build(args.output, args.shared, sid=args.sid, cache=args.cache)
        # The hash to pin as PACKAGED_MANIFEST_SHA256.
        print("%s  %s" % (_digest(args.output / "manifest.json"), args.output / "manifest.json"))
    elif args.command == "check":
        problems = check(args.shared, sid=args.sid, cache=args.cache)
        for problem in problems:
            print(problem, file=sys.stderr)
        if problems:
            return 1
        print("osteosarc makes Isovar's packaged bundle again from %s" % args.shared)
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
    sys.exit(main())
