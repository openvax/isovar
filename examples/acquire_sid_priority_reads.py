"""Acquire original regional RNA subsets for vaccine claims and selected fusions.

python -m examples.acquire_sid_priority_reads AUDIT --metadata-cache CACHE --cache CACHE
The AUDIT needs the frozen SV inventory and references. Whole libraries are never
downloaded. Vaccine claims are historical metadata, not a new treatment ranking.
"""

import argparse
from concurrent.futures import ProcessPoolExecutor, as_completed
from dataclasses import asdict
import multiprocessing
import os
from pathlib import Path
import shutil
import signal
import subprocess

from osteosarc import Cache, Dataset, RecoveryPolicy, extract_reads, read_receipt_files
from osteosarc.corpus import AUDIT_SNVS
from osteosarc.errors import OsteosarcError, RecordLimitError
from osteosarc.models import Variant

from examples.sid_sv_audit import acquisition, reconstruct, references
from examples.sid_sv_audit.inventory import digest, identity, load_inventory, read_json, write_json


FUSIONS = {"TPST1::CRCP", "GABBR1::SLC29A1", "PARD3B::CDKN2B-AS1;CDKN2B", "OTUD7A::FMN1"}
PRODUCT_PREFIXES = ("8275c5f21164", "58f81602ec2f", "40a3ceff825d", "6e9147690525")


def select_targets(variants, manifest):
    """Union vaccine claims, keep controls, and preserve distinct fusion geometries."""
    selected = {}
    for variant in variants:
        claims = dict(index_count=variant.vaccine_count, overlap=list(variant.vaccines),
                      source_variants=list(variant.annotations.get("source_vaccines", ())))
        vaccinated = bool((variant.vaccine_count or 0) > 0 or claims["overlap"] or claims["source_variants"])
        if vaccinated or variant.id in AUDIT_SNVS:
            selected[variant.id] = dict(variant=asdict(variant), vaccine_claims=claims,
                                       reasons=(["vaccine_membership"] if vaccinated else []) +
                                               (["RNA_control"] if variant.id in AUDIT_SNVS else []))
    geometries = {row["geometry_id"] for row in manifest["nominations"].values()
                  if row["original"].get("name") in FUSIONS and row["geometry_id"] is not None}
    return dict(small_variants=selected, geometries=sorted(geometries))


def acquire_variant_batch(directory, source, targets, cache, policy, selection_id, padding=150):
    """Partition oversized requests while keeping every target and original record."""
    names = sorted(targets)
    request = dict(source_identity=identity(source), targets=names, selection_id=selection_id,
                   padding=padding, policy=asdict(policy), timeout=600)
    path = Path(directory) / "small-variants" / source["id"] / (identity(request)[:24] + ".json")
    if path.exists():
        saved = read_json(path)
        if saved["request"] != request:
            raise ValueError("Focused acquisition request drift")
        acquisition.verify_upstream(saved)
        for key in ("bam", "index"):
            if key in saved and digest(saved[key]["path"]) != saved[key]["sha256"]:
                raise ValueError("Focused acquired input checksum mismatch")
        if saved["status"] == "partitioned":
            for group in saved["children"]:
                child = acquire_variant_batch(directory, source, {n: targets[n] for n in group["targets"]},
                                              cache, policy, selection_id, padding)
                if identity(child) != group["receipt_sha256"]:
                    raise ValueError("Focused partition child changed")
        return saved
    if shutil.disk_usage(cache.root).free < 8 * 1024 ** 3:
        raise OSError("Focused acquisition paused: less than 8 GiB free")
    result = dict(request=request, source_id=source["id"])
    header = acquisition.header(directory, source, cache)
    if header["status"] != "ready":
        result.update(status=header["status"])
    else:
        regions = [Variant(**targets[name]["variant"]).region(padding=padding) for name in names]
        try:
            subset = extract_reads(acquisition.source_file(source), regions, cache=cache,
                                   recovery=policy, timeout=600)
        except (OSError, ValueError, OsteosarcError, subprocess.SubprocessError) as error:
            if isinstance(error, (RecordLimitError, subprocess.TimeoutExpired, subprocess.CalledProcessError)) \
                    and len(names) > 1:
                children = []
                for group in (names[:len(names) // 2], names[len(names) // 2:]):
                    child = acquire_variant_batch(directory, source, {n: targets[n] for n in group},
                                                  cache, policy, selection_id, padding)
                    children.append(dict(targets=group, receipt_sha256=identity(child)))
                result.update(status="partitioned", children=children)
            else:
                result.update(status="acquisition_error")
            result["error"] = dict(type=type(error).__name__, message=str(error))
        else:
            result.update(status=subset.receipt.get("status", "acquired"),
                          receipt=acquisition.receipt_summary(subset.receipt),
                          upstream_receipt_files=[dict(path=str(p), sha256=h) for p, h in
                              read_receipt_files(subset.receipt_path, subset.receipt_path.parents[2]).items()],
                          bam=dict(path=str(subset.path), sha256=digest(subset.path)),
                          index=dict(path=str(subset.index_path), sha256=digest(subset.index_path)))
    write_json(path, result)
    return result


def run_source(directory, cache_root, source_id, selection_id, max_seconds):
    manifest, reference = load_inventory(directory), references.load(directory)
    selection = read_json(Path(directory) / "priority-selection.json")
    source = manifest["sources"][source_id]
    cache = Cache(cache_root)
    policy = RecoveryPolicy(max_rounds=4, max_intervals=20000, max_bases=20_000_000,
                            max_records=1_000_000, on_timeout="incomplete",
                            partner_batch_size=64, max_partner_queries=128, partner_timeout=30)
    ready = {name: row for name, row in selection["small_variants"].items()
             if row["variant"]["status"] == "ready" and row["variant"]["assembly"] == "GRCh38"}
    variants = acquire_variant_batch(directory, source, ready, cache, policy, selection_id,
                                     selection["variant_padding"])
    print(source_id[:12], "vaccine/control reads", variants["status"], flush=True)
    batches = {b["geometries"][0]: b for b in read_json(Path(directory) / "batches.json")}
    engine = reconstruct.implementation_identity()
    for gid in selection["geometries"]:
        acquired = acquisition.acquire_batch(directory, source, batches[gid], cache, policy, timeout=600,
                                              manifest=manifest, reference=reference, min_free_gib=8)
        print(source_id[:12], gid, "acquisition", acquired["status"], flush=True)
        result = reconstruct.run_batch(directory, source_id, batches[gid]["id"],
                                       reconstruct.PARAMETERS, max_seconds, engine)
        print(source_id[:12], gid, result["outcomes"], flush=True)


def run_selected_sources(directory, cache, sources, selection_id, seconds, workers):
    """Stop owned workers on interruption/failure; unfinished targets stay pending."""
    with ProcessPoolExecutor(max_workers=workers) as pool:
        try:
            futures = [pool.submit(run_source, directory, cache, sid, selection_id, seconds) for sid in sources]
            for future in as_completed(futures):
                future.result()
        except BaseException:
            children = multiprocessing.active_children()
            for child in children:
                if child.is_alive():
                    os.kill(child.pid, signal.SIGINT)
            for child in children:
                child.join(timeout=1)
                if child.is_alive():
                    child.terminate()
            raise


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("directory", type=Path)
    parser.add_argument("--metadata-cache", required=True)
    parser.add_argument("--cache", required=True)
    parser.add_argument("--workers", type=int, default=2)
    parser.add_argument("--max-seconds", type=int, default=3600)
    parser.add_argument("--variant-padding", type=int, default=150)
    args = parser.parse_args()
    if not 1 <= args.workers <= 4 or not 1 <= args.max_seconds <= 3600:
        parser.error("Use 1 through 4 workers and 1 through 3600 seconds per view")
    if args.variant_padding < 0:
        parser.error("Variant padding must be nonnegative")
    manifest, reference = load_inventory(args.directory), references.load(args.directory)
    data = Dataset.open(cache=Cache(args.metadata_cache, offline=True))
    selection = select_targets(data.variants(set="all"), manifest)
    selection.update(metadata_snapshot_id=data.manifest["id"], inventory_sha256=digest(args.directory / "inventory.json.gz"),
                     source_ids=sorted(sid for sid in manifest["sources"] if sid.startswith(PRODUCT_PREFIXES)),
                     variant_padding=args.variant_padding,
                     scope="selected_vaccine_claims_controls_and_fusions_four_platform_products")
    if len(selection["source_ids"]) != 4:
        raise ValueError("Expected all four pinned pilot products")
    write_json(args.directory / "priority-selection.json", selection)
    batches = [acquisition.make_batch(manifest, reference, [gid], 2000) for gid in sorted(manifest["geometries"])]
    write_json(args.directory / "batches.json", batches)
    Path(args.cache).mkdir(parents=True, exist_ok=True)
    run_selected_sources(args.directory, args.cache, selection["source_ids"], identity(selection),
                         args.max_seconds, args.workers)


if __name__ == "__main__":
    main()
