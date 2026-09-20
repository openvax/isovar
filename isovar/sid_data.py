"""Reproduce and use the small Sid read fixtures, without network during tests.

Regeneration requires ``isovar[data]`` and a named osteosarc snapshot. Only
indexed regions are acquired. The recipe selects complete original SAM records
by checksum; regional background reads never enter the distributable bundle.
"""

import argparse
from collections import Counter
import gzip
from functools import lru_cache
from hashlib import sha256
import json
import os
from pathlib import Path
import shutil


DATA = Path(__file__).with_name("data")
RECIPE = DATA / "sid-read-recipe.json.gz"
BUNDLE = DATA / "sid-reads.json.gz"


def read_json(path):
    with gzip.open(path, "rt") if str(path).endswith(".gz") else open(path) as handle:
        return json.load(handle)


def write_json(path, value):
    """Write deterministic gzip bytes, independent of path and clock."""
    with Path(path).open("xb") as raw:
        with gzip.GzipFile(filename="", mode="wb", fileobj=raw, mtime=0) as handle:
            handle.write(json.dumps(value, sort_keys=True, separators=(",", ":")).encode())


def sam_digest(line):
    return sha256(line.encode("ascii")).hexdigest()


def open_dataset(snapshot=None, cache=None, offline=False):
    """Open an explicitly created snapshot; never create one implicitly."""
    name = snapshot or os.environ.get("ISOVAR_SID_SNAPSHOT")
    if not name:
        raise ValueError("Specify --snapshot or ISOVAR_SID_SNAPSHOT (an osteosarc snapshot name)")
    return _open_dataset(name, cache or os.environ.get("ISOVAR_SID_CACHE"), offline)


@lru_cache(maxsize=4)
def _open_dataset(name, cache, offline):
    try:
        from osteosarc import Cache, Dataset
    except ImportError as error:
        raise RuntimeError("Regeneration requires Python >=3.10 and pip install 'isovar[data]'") from error
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
        asset = dataset.asset(url)
        path = dataset.download(asset)
        receipt = dict(url=url, sha256=sha256(path.read_bytes()).hexdigest(),
                       size=path.stat().st_size, snapshot_id=dataset.id, asset_id=asset.id)
    else:
        receipt = dataset.cache.fetch(url, max_bytes=256_000_000)
        path = dataset.cache.path(receipt)
    return path, receipt.to_dict() if hasattr(receipt, "to_dict") else receipt


def sam_regions(regions, assembly, reference_lengths=None):
    """Convert explicit 1-based inclusive SAM intervals to osteosarc Regions."""
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


def extract_regions(url, regions, assembly, output=None, *, dataset=None,
                    reference_lengths=None, timeout=600):
    """Delegate header/index/range acquisition to the snapshot's exact Asset."""
    dataset = dataset or open_dataset()
    asset = dataset.asset(url)
    subset = dataset.extract_reads(
        asset, sam_regions(regions, assembly, reference_lengths), timeout=timeout)
    if output is not None:
        output = Path(output)
        if output.exists() or Path(str(output) + ".bai").exists():
            raise FileExistsError(output)
        shutil.copyfile(subset.path, output)
        shutil.copyfile(subset.index_path, str(output) + ".bai")
    return subset


def asset_identity(asset):
    return {key: getattr(asset, key) for key in ("id", "key", "url", "size", "modified")}


def required_counts(recipe, source_id):
    """Repeated fixture use is shared; true source duplicates remain required."""
    counts = Counter()
    for fixture in recipe["fixtures"].values():
        if fixture["source"] == source_id:
            needed = Counter(fixture["records"])
            if fixture.get("record_references", False):
                needed = Counter(dict.fromkeys(needed, 1))
            counts |= needed
    return counts


def select_records(records, required):
    """Fail if even one required original record or duplicate is absent."""
    observed = Counter()
    selected = {}
    for record in records:
        line = record.to_string()
        checksum = sam_digest(line)
        if checksum in required:
            observed[checksum] += 1
            selected[checksum] = line
    missing = required - observed
    if missing:
        raise ValueError("Source is missing %d required SAM records: %s" % (
            sum(missing.values()), dict(missing)))
    return selected


def minimal_header(header, records):
    """Retain decoding identity, without shipping thousands of irrelevant PGs."""
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
    # Some source records themselves have RG tags absent from their header.
    # Preserve that source condition instead of inventing read-group metadata.
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


def generate(recipe, output, dataset, *, source_ids=None):
    """Regenerate only the exact test selections from indexed osteosarc data.

    Separate source files make an interrupted acquisition resumable. Every
    reused file is checked against the recipe and source snapshot identity.
    """
    from osteosarc import Region
    output = Path(output)
    output.mkdir(parents=True, exist_ok=True)
    chosen = sorted(recipe["sources"] if source_ids is None else source_ids)
    for identity in chosen:
        source = recipe["sources"][identity]
        asset = dataset.asset(source["asset"]["key"])
        if asset_identity(asset) != source["asset"]:
            raise ValueError("Snapshot asset changed: " + asset.key)
        required = required_counts(recipe, identity)
        path = output / (identity + ".json.gz")
        if path.exists():
            cached = read_json(path)
            validate_source(cached, required)
            if cached["asset"] != source["asset"] or cached["snapshot_id"] != dataset.id:
                raise ValueError("Acquisition snapshot/asset mismatch: " + identity)
            continue
        regions = [Region(**region) for region in source["regions"]]
        if not regions:
            raise ValueError("No bounded query regions: " + identity)
        subset = dataset.extract_reads(asset, regions)
        with subset.open() as handle:
            records = select_records(handle, required)
            header = minimal_header(handle.header.to_dict(), records.values())
        value = dict(asset=source["asset"], snapshot_id=dataset.id,
                     header=header, records=records,
                     duplicate_counts={k: n for k, n in required.items() if n > 1},
                     acquisition=subset.receipt)
        write_json(path, value)
        print(identity, len(records), "selected records", flush=True)
    return chosen


def validate_source(source, required):
    records = source["records"]
    if set(records) != set(required) or source["duplicate_counts"] != {
            k: n for k, n in required.items() if n > 1}:
        raise ValueError("Selected records differ from the test recipe")
    if any(sam_digest(line) != checksum for checksum, line in records.items()):
        raise ValueError("Original SAM checksum mismatch")


def pack(recipe, acquired, output):
    sources = {}
    for identity in sorted(recipe["sources"]):
        source = read_json(Path(acquired) / (identity + ".json.gz"))
        validate_source(source, required_counts(recipe, identity))
        if source["asset"] != recipe["sources"][identity]["asset"]:
            raise ValueError("Source asset differs from recipe")
        sources[identity] = source
    bundle = dict(schema_version=1, sources=sources, recipe_sha256=recipe_digest(recipe))
    write_json(output, bundle)
    return bundle


def verify(recipe=None, bundle=None):
    recipe = read_json(RECIPE) if recipe is None else recipe
    bundle = read_json(BUNDLE) if bundle is None else bundle
    if bundle["recipe_sha256"] != recipe_digest(recipe):
        raise ValueError("Bundle recipe checksum mismatch")
    if set(bundle["sources"]) != set(recipe["sources"]):
        raise ValueError("Bundle source set differs from recipe")
    for identity, source in bundle["sources"].items():
        validate_source(source, required_counts(recipe, identity))
        if source["asset"] != recipe["sources"][identity]["asset"]:
            raise ValueError("Bundle source identity differs from recipe")
    return dict(fixtures=len(recipe["fixtures"]), sources=len(bundle["sources"]),
                unique_records=sum(len(s["records"]) for s in bundle["sources"].values()))


def recipe_digest(recipe):
    return sha256(json.dumps(recipe, sort_keys=True, separators=(",", ":")).encode()).hexdigest()


def export_fixture(name, output, *, recipe=None, bundle=None):
    """Export original SAM records in fixture order, including repeated records."""
    import pysam
    recipe = read_json(RECIPE) if recipe is None else recipe
    bundle = read_json(BUNDLE) if bundle is None else bundle
    verify(recipe, bundle)
    fixture = recipe["fixtures"][name]
    source = bundle["sources"][fixture["source"]]
    with Path(output).open("x") as handle:
        handle.write(str(pysam.AlignmentHeader.from_dict(source["header"])))
        for checksum in fixture["records"]:
            handle.write(source["records"][checksum] + "\n")


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("command", choices=("list", "verify", "generate", "pack", "export"))
    parser.add_argument("--recipe", type=Path, default=RECIPE)
    parser.add_argument("--bundle", type=Path, default=BUNDLE)
    parser.add_argument("--snapshot")
    parser.add_argument("--cache")
    parser.add_argument("--offline", action="store_true")
    parser.add_argument("--source", action="append")
    parser.add_argument("--acquired", type=Path)
    parser.add_argument("--output", type=Path)
    parser.add_argument("--fixture")
    args = parser.parse_args(argv)
    recipe = read_json(args.recipe)
    if args.command == "list":
        for name, fixture in recipe["fixtures"].items():
            print(name, len(fixture["records"]))
    elif args.command == "verify":
        print(json.dumps(verify(recipe, read_json(args.bundle)), sort_keys=True))
    elif args.command == "generate":
        if args.acquired is None:
            parser.error("generate requires --acquired")
        generate(recipe, args.acquired, open_dataset(args.snapshot, args.cache, args.offline),
                 source_ids=args.source)
    elif args.command == "pack":
        if args.acquired is None or args.output is None:
            parser.error("pack requires --acquired and --output")
        pack(recipe, args.acquired, args.output)
    else:
        if args.fixture is None or args.output is None:
            parser.error("export requires --fixture and --output")
        export_fixture(args.fixture, args.output, recipe=recipe, bundle=read_json(args.bundle))


if __name__ == "__main__":
    main()
