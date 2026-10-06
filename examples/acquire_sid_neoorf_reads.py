"""Acquire additional Sid indel, splice-site and documented fusion RNA inputs.

python -m examples.acquire_sid_neoorf_reads PRIOR_AUDIT ADDITIONS \
    --metadata-cache CACHE --cache ADDITIONS/cache
The prior audit is immutable. These inputs nominate possible non-SNV ORFs; they
do not establish altered splicing, protein expression or DNA-SV origin.
"""

import argparse
from copy import deepcopy
from pathlib import Path
import shutil

from osteosarc import Cache, Dataset

from examples import acquire_sid_priority_reads as priority
from examples.sid_sv_audit import acquisition, references
from examples.sid_sv_audit.inventory import digest, geometry, identity, load_inventory, read_json, write_json


FIXTURE = Path(__file__).resolve().parents[1] / "tests/data/fusions/coding-corpus/ATP5MG--KMT2A.input.json.gz"


def supplemental_inventory(manifest, fixture, fixture_sha256, parent_sha256):
    """Derive retained sides from the pinned original RNA's interbase breakpoints."""
    fusion = fixture["fusion"]
    if fusion["reference_name"] != "GRCh38" or fusion["event_id"] != "ATP5MG--KMT2A":
        raise ValueError("Unexpected supplemental fusion or assembly")
    ends = []
    for role in ("donor", "acceptor"):
        end = fusion[role]
        if end["strand"] not in ("+", "-"):
            raise ValueError("Unresolved supplemental orientation")
        left = (role == "donor") == (end["strand"] == "+")
        ends.append(dict(contig=end["contig"], position=end["position"],
                         retained_side="left" if left else "right"))
    ends = geometry(ends)
    gid = "adj-" + identity(ends)[:24]
    name = "fixture:ATP5MG--KMT2A"
    result = deepcopy(manifest)
    if name in result["nominations"]:
        raise ValueError("Supplemental nomination already exists")
    group = result["geometries"].setdefault(gid, dict(breakends=ends, nominations=[]))
    if group["breakends"] != ends:
        raise ValueError("Supplemental geometry collision")
    group["nominations"].append(name)
    result["nominations"][name] = dict(catalogue="original_RNA_fixture", geometry_id=gid,
        original=dict(name="ATP5MG::KMT2A", fusion=fusion, fixture_sha256=fixture_sha256,
                      note=fixture["note"]))
    result["supplemental_inputs"] = [dict(path="supplemental-inputs/ATP5MG--KMT2A.input.json.gz",
                                        sha256=fixture_sha256)]
    result["parent_inventory"] = dict(path="parent-inventory.json.gz", sha256=parent_sha256)
    return result


def prepare(prior, directory, metadata_cache):
    """Freeze only additions, retaining parent and supplemental source pins."""
    directory.mkdir(parents=True, exist_ok=True)
    manifest = supplemental_inventory(load_inventory(prior), read_json(FIXTURE), digest(FIXTURE),
                                      digest(prior / "inventory.json.gz"))
    write_json(directory / "inventory.json.gz", manifest)
    write_json(directory / "inventory-pin.json", dict(sha256=digest(directory / "inventory.json.gz")))
    for src, dest in [(prior / "inventory.json.gz", directory / "parent-inventory.json.gz"),
                      (FIXTURE, directory / manifest["supplemental_inputs"][0]["path"])]:
        dest.parent.mkdir(parents=True, exist_ok=True)
        if dest.exists():
            if digest(src) != digest(dest):
                raise ValueError("Derived inventory source changed")
        else:
            shutil.copyfile(src, dest)
    data = Dataset.open(cache=Cache(metadata_cache, offline=True))
    previous = read_json(prior / "priority-selection.json")
    if data.manifest["id"] != previous["metadata_snapshot_id"]:
        raise ValueError("Use the same metadata snapshot as the parent audit")
    selection = priority.select_targets(data.variants(set="all"), manifest, include_neo_orfs=True)
    selection["small_variants"] = {k: v for k, v in selection["small_variants"].items()
                                   if k not in previous["small_variants"]}
    selection["geometries"] = sorted(set(selection["geometries"]) - set(previous["geometries"]))
    selection.update(metadata_snapshot_id=data.manifest["id"],
                     inventory_sha256=digest(directory / "inventory.json.gz"),
                     source_ids=previous["source_ids"], variant_padding=previous["variant_padding"],
                     parent_selection_sha256=digest(prior / "priority-selection.json"),
                     scope="additional_indel_splice_and_documented_fusion_inputs_four_products")
    write_json(directory / "priority-selection.json", selection)
    write_json(directory / "parent-priority-selection.json", previous)
    if not (directory / "references-pin.json").exists():
        reference = references.build(directory)
    else:
        reference = references.load(directory)
    batches = [acquisition.make_batch(manifest, reference, [gid], 2000) for gid in selection["geometries"]]
    write_json(directory / "batches.json", batches)
    return selection


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("prior", type=Path)
    parser.add_argument("directory", type=Path)
    parser.add_argument("--metadata-cache", required=True)
    parser.add_argument("--cache", required=True)
    parser.add_argument("--workers", type=int, default=2)
    parser.add_argument("--max-seconds", type=int, default=3600)
    args = parser.parse_args()
    if not 1 <= args.workers <= 4 or not 1 <= args.max_seconds <= 3600:
        parser.error("Use 1 through 4 workers and 1 through 3600 seconds per view")
    if args.prior.resolve() == args.directory.resolve():
        parser.error("Use a new directory; the prior audit is immutable")
    selection = prepare(args.prior, args.directory, args.metadata_cache)
    Path(args.cache).mkdir(parents=True, exist_ok=True)
    priority.run_selected_sources(args.directory, args.cache, selection["source_ids"],
                                  identity(selection), args.max_seconds, args.workers)


if __name__ == "__main__":
    main()
