"""Independently verify all ten intended dense-pilot view outcomes and witnesses.

python -m examples.check_sid_sv_dense AUDIT OUTPUT
"""

import argparse
from collections import Counter
from hashlib import sha256
from pathlib import Path
import subprocess
import tempfile

import pysam

from examples.check_sid_sv_pilot import check_candidate
from examples.rerun_sid_sv_dense import SOURCE, TARGETS
from examples.sid_sv_audit import acquisition, reconstruct, references
from examples.sid_sv_audit.inventory import digest, identity, load_inventory, read_json, write_json


def verify(directory, output):
    manifest, reference = load_inventory(directory), references.load(directory)
    batches = read_json(directory / "batches.json")
    by_geometry = {gid: batch for batch in batches for gid in batch["geometries"]}
    rows, checks, outcomes, pins = [], [], Counter(), {}
    engine = None
    for gid, name in TARGETS.items():
        batch = by_geometry[gid]
        leaves, receipts = acquisition.acquisition_tree(directory, SOURCE, batch["id"])
        assert len(leaves) == 1
        receipt = next(iter(leaves.values()))
        assert receipt["status"] in acquisition.USABLE_STATUSES
        assert receipt["request"]["source_identity"] == identity(manifest["sources"][SOURCE])
        for key in ("bam", "index"):
            assert digest(receipt[key]["path"]) == receipt[key]["sha256"]
        for bid in receipts:
            path = directory / "sources" / SOURCE / (bid + ".json")
            pins[str(path.relative_to(directory))] = digest(path)
        bam_path = receipt["bam"]["path"]
        native_total = int(subprocess.run(["samtools", "view", "-c", bam_path],
                                         capture_output=True, check=True).stdout)
        with tempfile.TemporaryDirectory() as temporary:
            bed = Path(temporary) / "regions.bed"
            bed.write_text("".join("%s\t%d\t%d\n" % tuple(region) for region in batch["regions"]))
            native_regional = int(subprocess.run(["samtools", "view", "-c", "-M", "-L", str(bed), bam_path],
                                                capture_output=True, check=True).stdout)
        views, candidates = {}, []
        for view in ("forward", "reverse"):
            path = reconstruct.result_path(directory, SOURCE, gid, view)
            result = read_json(path)
            request = result["request"]
            assert request["parameters"]["dense_support"]
            assert request["acquisition_sha256"] == identity(receipt)
            assert request["inventory_sha256"] == digest(directory / "inventory.json.gz")
            assert request["references_sha256"] == digest(directory / "references.json.gz")
            assert request["source_id"] == SOURCE and request["geometry_id"] == gid
            assert request["orientation"] == view
            if engine is not None:
                assert engine == request["engine_id"]
            engine = request["engine_id"]
            # A timeout or missing result must not masquerade as a full scan.
            assert "reconstruction" in result
            support = result["support_acquisition"]
            assert support["records_complete"] and support["records_scanned"] == native_total
            details = directory / result["reconstruction"]["path"]
            assert digest(details) == result["reconstruction"]["sha256"]
            pins.update({str(p.relative_to(directory)): digest(p) for p in (path, details)})
            views[view] = {key: result[key] for key in (
                "status", "limitations", "records", "discovery", "support_acquisition", "input_limitations")}
            views[view]["candidates"] = len(result["candidates"])
            outcomes[result["status"]] += 1
            candidates.extend(dict(c, source_id=SOURCE, geometry_id=gid, orientation=view)
                              for c in result["candidates"])
        if candidates:
            candidate = max(candidates, key=lambda c: c["rna_support"]["fragments"])
            details = read_json(directory / "reconstructions" / SOURCE / (
                gid + "." + candidate["orientation"] + ".json.gz"))
            check = check_candidate(candidate, details, reference)
            expected = set(check["original_record_hashes"].values())
            names = {details["original_records"][key].split("\t", 1)[0]
                     for key in check["original_record_hashes"]}
            observed = set()
            with pysam.AlignmentFile(bam_path) as bam:
                for read in bam:
                    if read.query_name in names:
                        observed.add(sha256(read.to_string().encode()).hexdigest())
            assert expected <= observed
            check["original_bam_records_verified"] = len(expected)
            checks.append(check)
        rows.append(dict(geometry_id=gid, target=name, nominations=manifest["geometries"][gid]["nominations"],
                         acquisition_status=receipt["status"], regional_alignment_records=native_regional,
                         original_input_records=native_total, bam_sha256=receipt["bam"]["sha256"], views=views))
    ledger = dict(scope="five_dense_10x_pilot_geometries", source_id=SOURCE,
                  intended_targets=len(TARGETS), intended_views=2 * len(TARGETS),
                  accounted_views=sum(outcomes.values()), outcome_counts=dict(outcomes),
                  inventory_sha256=digest(directory / "inventory.json.gz"),
                  references_sha256=digest(directory / "references.json.gz"), engine_id=engine,
                  rows=rows, independent_checks=checks, pins=pins)
    assert ledger["accounted_views"] == ledger["intended_views"]
    write_json(output / "dense-pilot.json.gz", ledger)
    print(dict(outcomes), "independently checked candidates:", len(checks))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("directory", type=Path)
    parser.add_argument("output", type=Path)
    args = parser.parse_args()
    verify(args.directory, args.output)


if __name__ == "__main__":
    main()
