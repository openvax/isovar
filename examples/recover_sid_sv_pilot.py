"""Recover every original acquisition-failed Sid pilot pair with pinned inputs.

Run acquisition and reconstruction separately or together. The selected scope
always includes every failed pair, even when --pair runs only one of them.
"""

import argparse
from collections import Counter
from dataclasses import asdict
from importlib.metadata import version
import json
import os
from pathlib import Path
import shutil
import tempfile

from osteosarc import Cache, RecoveryPolicy
from packaging.version import Version

import isovar
from examples.sid_sv_audit import acquisition, reconstruct, references
from examples.sid_sv_audit.inventory import digest, identity, load_inventory, read_json, write_json


def failed_pairs(pilot):
    """Retain both original failed views and their exact source/geometry identity."""
    pairs = {}
    for row in pilot["outcomes"]:
        if row["status"] == "not_assessable":
            key = row["source_id"] + "/" + row["geometry_id"]
            pairs.setdefault(key, []).append(row)
    if not pairs or any(sorted(r["orientation"] for r in rows) != ["forward", "reverse"]
                        for rows in pairs.values()):
        raise ValueError("Pilot must contain both failed orientations for every selected pair")
    return [dict(pair=key, source_id=rows[0]["source_id"], geometry_id=rows[0]["geometry_id"],
                 original_outcomes=sorted(rows, key=lambda r: r["orientation"]))
            for key, rows in sorted(pairs.items())]


def prepare(directory, base_audit, pilot_path):
    """Freeze a selected denominator; an existing directory must match exactly."""
    directory, base_audit, pilot_path = map(Path, (directory, base_audit, pilot_path))
    manifest, reference = load_inventory(base_audit), references.load(base_audit)
    pilot = read_json(pilot_path)
    if (pilot["inventory_sha256"] != digest(base_audit / "inventory.json.gz")
            or pilot["references_sha256"] != digest(base_audit / "references.json.gz")):
        raise ValueError("Pilot and frozen inventory/reference inputs differ")
    pairs = failed_pairs(pilot)
    gids = sorted({row["geometry_id"] for row in pairs})
    if any(row["source_id"] not in manifest["sources"] for row in pairs):
        raise ValueError("Pilot source is absent from the inventory")
    batches = [acquisition.make_batch(manifest, reference, [gid], 2000) for gid in gids]
    selection = dict(scope="original_acquisition_failed_pilot_pairs", pairs=pairs,
                     original_pilot_sha256=digest(pilot_path),
                     original_pilot_views=pilot["expected_geometry_views"],
                     inventory_sha256=pilot["inventory_sha256"],
                     references_sha256=pilot["references_sha256"])
    files = ("inventory.json.gz", "inventory-pin.json", "references.json.gz", "references-pin.json")
    if directory.exists():
        if read_json(directory / "recovery-selection.json") != selection:
            raise ValueError("Recovery selection drift")
        if read_json(directory / "batches.json") != batches:
            raise ValueError("Recovery batch drift")
        for name in files:
            if digest(directory / name) != digest(base_audit / name):
                raise ValueError("Recovery frozen input changed: " + name)
    else:
        directory.parent.mkdir(parents=True, exist_ok=True)
        with tempfile.TemporaryDirectory(dir=directory.parent, prefix=".sid-recovery-") as temporary:
            work = Path(temporary)
            for name in files:
                shutil.copyfile(base_audit / name, work / name)
            write_json(work / "batches.json", batches)
            write_json(work / "recovery-selection.json", selection)
            work.rename(directory)
    return manifest, reference, selection, {b["geometries"][0]: b for b in batches}


def clone_cache(source, destination):
    """Copy metadata separately and share immutable BAM/index bytes by hard link."""
    source, destination = Path(source).resolve(), Path(destination).resolve()
    if source == destination or source in destination.parents or destination in source.parents:
        raise ValueError("Reuse cache must be a separate directory")
    if destination.exists():
        return

    def copy(path, target):
        if Path(path).suffix in (".bam", ".bai"):
            os.link(path, target)
            return target
        return shutil.copy2(path, target)

    destination.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(dir=destination.parent, prefix=".seed-cache-") as temporary:
        work = Path(temporary) / "cache"
        shutil.copytree(source, work, copy_function=copy)
        work.rename(destination)


def coverage(directory, selection, batches):
    """Account for every selected view; absent output is pending, never negative."""
    directory = Path(directory)
    rows = []
    for pair in selection["pairs"]:
        sid, gid = pair["source_id"], pair["geometry_id"]
        acquired_path = directory / "sources" / sid / (batches[gid]["id"] + ".json")
        acquired = read_json(acquired_path) if acquired_path.exists() else None
        for orientation in ("forward", "reverse"):
            path = reconstruct.result_path(directory, sid, gid, orientation)
            row = dict(source_id=sid, geometry_id=gid, orientation=orientation,
                       acquisition_status=acquired["status"] if acquired else "pending",
                       status="pending_reconstruction", reconstructed=False,
                       original_outcomes=pair["original_outcomes"])
            if acquired:
                row["acquisition_sha256"] = digest(acquired_path)
                row["acquired_records"] = acquired.get("receipt", {}).get("records")
                row["acquisition_scope"] = acquired.get("receipt", {}).get("scope")
            if path.exists():
                result = read_json(path)
                request = result["request"]
                if (request["source_id"], request["geometry_id"], request["orientation"]) != (sid, gid, orientation):
                    raise ValueError("Recovery result is filed under the wrong identity")
                row.update(status=result["status"], result_sha256=digest(path),
                           result_path=str(path.resolve()), reason=result.get("reason"),
                           candidate_sequences=len(result["candidates"]) if "candidates" in result else None,
                           input_limitations=result["input_limitations"],
                           discovery=result.get("discovery"), support_acquisition=result.get("support_acquisition"))
                if "reconstruction" in result:
                    details = (directory / result["reconstruction"]["path"]).resolve()
                    if (not details.is_relative_to(directory.resolve())
                            or digest(details) != result["reconstruction"]["sha256"]):
                        raise ValueError("Missing or changed completed reconstruction")
                    row.update(reconstructed=True, reconstruction_sha256=result["reconstruction"]["sha256"])
            rows.append(row)
    return dict(scope=selection["scope"], selection_sha256=identity(selection), outcomes=rows,
                expected_pairs=len(selection["pairs"]), expected_views=2 * len(selection["pairs"]),
                acquisition_counts=dict(Counter(row["acquisition_status"] for row in rows)),
                outcome_counts=dict(Counter(row["status"] for row in rows)),
                all_inputs_acquired=all(row["acquisition_status"] in acquisition.USABLE_STATUSES for row in rows),
                all_views_reconstructed=all(row["reconstructed"] for row in rows), full_catalogue_complete=False)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("directory", type=Path)
    parser.add_argument("--base-audit", required=True, type=Path)
    parser.add_argument("--pilot-coverage", required=True, type=Path)
    parser.add_argument("--cache", required=True, type=Path)
    parser.add_argument("--reuse-cache", type=Path)
    parser.add_argument("--stage", choices=("acquire", "reconstruct", "all"), default="all")
    parser.add_argument("--pair", action="append", help="Exact source-id/geometry-id; denominator stays complete")
    parser.add_argument("--max-acquisition-records", type=int, default=6_500_000)
    parser.add_argument("--seed-timeout", type=int, default=1200)
    parser.add_argument("--max-seconds", type=int, default=3600)
    parser.add_argument("--min-free-gib", type=int, default=8)
    args = parser.parse_args()
    if min(args.max_acquisition_records, args.seed_timeout, args.max_seconds) < 1 or args.min_free_gib < 0:
        parser.error("Record/time budgets must be positive; the free-space guard cannot be negative")
    if Version(version("osteosarc")) < Version("0.15.5"):
        parser.error("Osteosarc 0.15.5 or newer is required for snapshot source-date checks")
    manifest, reference, selection, batches = prepare(args.directory, args.base_audit, args.pilot_coverage)
    requested = set(args.pair or [row["pair"] for row in selection["pairs"]])
    if requested - {row["pair"] for row in selection["pairs"]}:
        parser.error("Unknown recovery pair")
    policy = RecoveryPolicy(max_rounds=4, max_intervals=20000, max_bases=20_000_000,
                            max_records=args.max_acquisition_records, on_timeout="incomplete",
                            partner_batch_size=64, max_partner_queries=128, partner_timeout=30)
    engine = reconstruct.implementation_identity()
    run = dict(selection_sha256=identity(selection), acquisition_policy=asdict(policy),
               seed_timeout=args.seed_timeout, max_seconds=args.max_seconds,
               parameters=reconstruct.PARAMETERS, engine_id=engine, isovar_version=isovar.__version__,
               osteosarc_version=version("osteosarc"), code_sha256=digest(__file__))
    run_path = args.directory / "recovery-run.json"
    if run_path.exists() and read_json(run_path) != run:
        raise ValueError("Recovery request drift; use a new attempt directory")
    write_json(run_path, run)
    if args.reuse_cache:
        clone_cache(args.reuse_cache, args.cache)
    args.cache.mkdir(parents=True, exist_ok=True)
    cache = Cache(args.cache)
    selected = [row for row in selection["pairs"] if row["pair"] in requested]
    for stage in ("acquire", "reconstruct"):
        if args.stage not in (stage, "all"):
            continue
        for row in selected:
            sid, gid = row["source_id"], row["geometry_id"]
            print(json.dumps(dict(stage=stage, pair=row["pair"], status="starting")), flush=True)
            if stage == "acquire":
                result = acquisition.acquire_batch(args.directory, manifest["sources"][sid], batches[gid],
                    cache, policy, timeout=args.seed_timeout, manifest=manifest,
                    reference=reference, min_free_gib=args.min_free_gib)
                print(json.dumps(dict(stage=stage, pair=row["pair"], status=result["status"],
                                      records=result.get("receipt", {}).get("records"))), flush=True)
            else:
                result = reconstruct.run_batch(args.directory, sid, batches[gid]["id"],
                                               reconstruct.PARAMETERS, args.max_seconds, engine)
                print(json.dumps(dict(stage=stage, pair=row["pair"], result=result)), flush=True)
    observed = coverage(args.directory, selection, batches)
    path = args.directory / ("recovery-coverage." + identity(observed) + ".json.gz")
    write_json(path, observed)
    print(json.dumps(dict(stage="coverage", path=str(path), expected_views=observed["expected_views"],
                          acquisition_counts=observed["acquisition_counts"],
                          outcome_counts=observed["outcome_counts"])), flush=True)


if __name__ == "__main__":
    main()
