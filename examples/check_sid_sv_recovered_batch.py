"""Verify and publish the recovered 13-target Sid 10x pilot batch offline.

python -m examples.check_sid_sv_recovered_batch ORIGINAL RETRY DENSE CONFIRMED OUTPUT
"""

import argparse
from collections import Counter
from hashlib import sha256
from pathlib import Path
import subprocess
import tempfile

import pysam

from examples.check_sid_sv_pilot import check_candidate
from examples.sid_sv_audit import acquisition, reconstruct, references
from examples.sid_sv_audit.inventory import digest, identity, load_inventory, read_json, write_json

SOURCE = "6e9147690525529f9591602afd1191f9862f4df16117a865c6726654456c9368"
BATCH = "8f04d8507513062972345a77"


def verify(original, retry, dense, confirmed, output):
    manifest, reference = load_inventory(confirmed), references.load(confirmed)
    baseline = original / "sources" / SOURCE / (BATCH + ".json")
    failure = read_json(baseline)
    assert failure["status"] == "acquisition_error"
    intended = failure["request"]["batch"]["geometries"]
    assert len(intended) == len(set(intended)) == 13
    acquired = {}
    for path in (confirmed / "sources" / SOURCE).glob("*.json"):
        if path.name == "header.json":
            continue
        value = read_json(path)
        acquisition.verify_upstream(value)
        assert value["request"]["source_identity"] == identity(manifest["sources"][SOURCE])
        for gid in value["request"]["batch"]["geometries"]:
            if gid in intended:
                assert gid not in acquired
                acquired[gid] = (path, value)
    assert set(acquired) == set(intended)
    rows, checks, pins = [], [], {str(baseline): digest(baseline)}
    outcomes = Counter()
    verified_bams = set()
    for gid in intended:
        path, receipt = acquired[gid]
        assert receipt["status"] == "acquired"
        pins[str(path)] = digest(path)
        bam_path = receipt["bam"]["path"]
        if bam_path not in verified_bams:
            assert digest(bam_path) == receipt["bam"]["sha256"]
            assert digest(receipt["index"]["path"]) == receipt["index"]["sha256"]
            verified_bams.add(bam_path)
        batch = acquisition.make_batch(manifest, reference, [gid], failure["request"]["batch"]["window"])
        with tempfile.TemporaryDirectory() as temporary:
            bed = Path(temporary) / "regions.bed"
            bed.write_text("".join("%s\t%d\t%d\n" % tuple(region) for region in batch["regions"]))
            command = ["samtools", "view", "-c", "-M", "-L", str(bed), bam_path]
            native_count = int(subprocess.run(command, capture_output=True, check=True).stdout)
        views, candidates = {}, []
        for view in ("forward", "reverse"):
            result_path = reconstruct.result_path(confirmed, SOURCE, gid, view)
            result = read_json(result_path)
            before = read_json(reconstruct.result_path(original, SOURCE, gid, view))
            assert before["status"] == "not_assessable"
            assert result["request"]["acquisition_sha256"] == identity(receipt)
            first = reconstruct.result_path(retry, SOURCE, gid, view)
            if not first.exists() or read_json(first)["status"] == "not_assessable":
                first = reconstruct.result_path(dense, SOURCE, gid, view)
            earlier = read_json(first)
            assert {k: v for k, v in earlier["request"].items() if k != "engine_id"} == {
                k: v for k, v in result["request"].items() if k != "engine_id"}
            assert {k: v for k, v in earlier.items() if k != "request"} == {
                k: v for k, v in result.items() if k != "request"}
            details = confirmed / result["reconstruction"]["path"]
            assert digest(details) == result["reconstruction"]["sha256"]
            pins.update({str(p): digest(p) for p in (result_path, first, details)})
            views[view] = {key: result[key] for key in ("status", "records", "limitations")}
            views[view]["candidates"] = len(result["candidates"])
            views[view].update(engine_id=result["request"]["engine_id"],
                               first_engine_id=earlier["request"]["engine_id"],
                               reconstruction_sha256=result["reconstruction"]["sha256"])
            outcomes[result["status"]] += 1
            candidates.extend(dict(c, source_id=SOURCE, geometry_id=gid, orientation=view)
                              for c in result["candidates"])
        if candidates:
            candidate = max(candidates, key=lambda c: c["rna_support"]["fragments"])
            details = read_json(confirmed / "reconstructions" / SOURCE / (gid + "." + candidate["orientation"] + ".json.gz"))
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
        rows.append(dict(geometry_id=gid, nominations=manifest["geometries"][gid]["nominations"],
                         original_status="not_assessable", acquisition_status=receipt["status"],
                         regional_alignment_records=native_count, views=views))
    result = dict(scope="original_failed_13_target_10x_batch", source_id=SOURCE, batch_id=BATCH,
                  intended_targets=13, intended_views=26, accounted_views=sum(outcomes.values()),
                  outcome_counts=dict(outcomes), offline_repeat_identical=True,
                  inventory_sha256=digest(confirmed / "inventory.json.gz"),
                  references_sha256=digest(confirmed / "references.json.gz"),
                  rows=rows, independent_checks=checks, pins=pins)
    assert result["accounted_views"] == 26
    write_json(output / "recovered-batch.json.gz", result)
    print(dict(outcomes), "independent candidates:", len(checks))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ("original", "retry", "dense", "confirmed", "output"):
        parser.add_argument(name, type=Path)
    args = parser.parse_args()
    verify(args.original, args.retry, args.dense, args.confirmed, args.output)


if __name__ == "__main__":
    main()
