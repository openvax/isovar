"""Acquire and replay the five dense Sid 10x pilot geometries in a fresh audit.

Preserves the complete catalogue inventory and batch plan. Only the selected
acquisitions run; the remaining catalogue stays pending in the coverage report.
"""

import argparse
from pathlib import Path

from osteosarc import Cache, RecoveryPolicy

from examples.sid_sv_audit import acquisition, reconstruct, references
from examples.sid_sv_audit.inventory import load_inventory, write_json


SOURCE = "6e9147690525529f9591602afd1191f9862f4df16117a865c6726654456c9368"
TARGETS = {
    "adj-02ce630f08dfbac9b6c4cde5": "::CNN2",
    "adj-041dbfbb2239d183a5982838": "::IGKC",
    "adj-0742149ff01fe14252c0cd26": "UBC::RPS27A",
    "adj-06dc27c7e4a9fbf7768d2c26": "FTH1::FTL",
    "adj-01b50c22a55d873dce244b71": "TMSB10::TMSB4X",
}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("directory", type=Path)
    parser.add_argument("--cache", required=True)
    parser.add_argument("--target", choices=TARGETS, action="append")
    parser.add_argument("--max-seconds", type=int, default=3600)
    args = parser.parse_args()
    if not 1 <= args.max_seconds <= 3600:
        parser.error("Per-view time limit must be 1 through 3600 seconds")
    manifest, reference = load_inventory(args.directory), references.load(args.directory)
    targets = args.target or list(TARGETS)
    # Single-target batches for this selected pilot; all other geometries remain
    # planned. This never shrinks the catalogue's intended denominator.
    batches = [acquisition.make_batch(manifest, reference, [gid], 2000)
               for gid in sorted(manifest["geometries"])]
    write_json(args.directory / "batches.json", batches)
    by_geometry = {b["geometries"][0]: b for b in batches}
    cache = Cache(args.cache)
    cache.root.mkdir(parents=True, exist_ok=True)
    policy = RecoveryPolicy(max_rounds=4, max_intervals=20000, max_bases=20_000_000,
                            max_records=6_500_000, on_timeout="incomplete",
                            partner_batch_size=64, max_partner_queries=128, partner_timeout=30)
    engine = reconstruct.implementation_identity()
    source = manifest["sources"][SOURCE]
    for gid in targets:
        print(TARGETS[gid], "acquiring", flush=True)
        receipt = acquisition.acquire_batch(args.directory, source, by_geometry[gid], cache, policy,
                                            timeout=600, manifest=manifest, reference=reference, min_free_gib=8)
        print(TARGETS[gid], "acquisition", receipt["status"], flush=True)
        result = reconstruct.run_batch(args.directory, SOURCE, by_geometry[gid]["id"],
                                       reconstruct.PARAMETERS, args.max_seconds, engine)
        print(TARGETS[gid], "reconstruction", result, flush=True)


if __name__ == "__main__":
    main()
