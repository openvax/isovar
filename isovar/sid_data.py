"""Reproduce and use the small Sid read fixtures, without network during tests.

Regeneration requires a named osteosarc snapshot. Only
indexed regions are acquired. The recipe selects complete original SAM records
by checksum; regional background reads never enter the distributable bundle.
"""

from osteosarc import legacy_fixtures as shared
from osteosarc.legacy_fixtures import (
    read_json as read_json,
    write_json as write_json,
    sam_digest as sam_digest,
    sam_regions as sam_regions,
    asset_identity as asset_identity,
    required_counts as required_counts,
    select_records as select_records,
    minimal_header as minimal_header,
    generate as generate,
    validate_source as validate_source,
    pack as pack,
    recipe_digest as recipe_digest,
)

import argparse
from functools import lru_cache
from hashlib import sha256
import json
import os
from pathlib import Path
import shutil


DATA = Path(__file__).with_name("data")
RECIPE = DATA / "sid-read-recipe.json.gz"
BUNDLE = DATA / "sid-reads.json.gz"


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


def verify(recipe=None, bundle=None):
    return shared.verify(read_json(RECIPE) if recipe is None else recipe,
                         read_json(BUNDLE) if bundle is None else bundle)


def export_fixture(name, output, *, recipe=None, bundle=None):
    return shared.export_fixture(name, output, recipe=read_json(RECIPE) if recipe is None else recipe,
                                 bundle=read_json(BUNDLE) if bundle is None else bundle)


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
    parser.add_argument("--panel-recipe", type=Path, help="Shared Osteosarc v1 panel recipe")
    parser.add_argument("--panel-source", action="append", default=[], metavar="ID=LOCAL_BAM")
    args = parser.parse_args(argv)
    if args.panel_recipe:
        if args.command != "generate" or args.output is None:
            parser.error("--panel-recipe requires generate and --output")
        from osteosarc import Cache
        from osteosarc.bundles import generate_panel
        generate_panel(args.panel_recipe, args.output, sources=args.panel_source,
                       cache=Cache(args.cache, offline=args.offline))
        return

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
