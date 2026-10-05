"""Verify the focused Sid subset denominator, original BAM pins and RNA witnesses."""

import argparse
from collections import Counter
import csv
from hashlib import sha256
from pathlib import Path
import subprocess

import pysam

from examples.check_sid_sv_pilot import check_candidate
from examples.sid_sv_audit import acquisition, reconstruct, references
from examples.sid_sv_audit.inventory import digest, identity, load_inventory, read_json, write_json


def verify(directory, output):
    manifest, reference = load_inventory(directory), references.load(directory)
    selection = read_json(directory / "priority-selection.json")
    selection_id = identity(selection)
    ready = {n for n, r in selection["small_variants"].items()
             if r["variant"]["status"] == "ready" and r["variant"]["assembly"] == "GRCh38"}
    if ready != set(selection["small_variants"]):
        raise ValueError("Selected unresolved/reference-mismatched alleles need explicit outcomes")
    rows, checks, pins, candidates = [], [], {}, []
    for filename in ("inventory.json.gz", "inventory-pin.json", "references.json.gz", "references-pin.json",
                     "priority-selection.json", "batches.json"):
        pins[filename] = digest(directory / filename)

    def pin(path):
        path = Path(path)
        pins[str(path.relative_to(directory))] = digest(path)

    def verify_input(receipt):
        acquisition.verify_upstream(receipt)
        for asset in receipt.get("upstream_receipt_files", []):
            pin(asset["path"])
        for key in ("bam", "index"):
            if key in receipt:
                assert digest(receipt[key]["path"]) == receipt[key]["sha256"]
                pin(receipt[key]["path"])

    for sid in selection["source_ids"]:
        source = manifest["sources"][sid]
        pin(directory / "sources" / sid / "header.json")
        files = list((directory / "small-variants" / sid).glob("*.json"))
        saved = {p.name: (p, read_json(p)) for p in files}
        roots = [(p, r) for p, r in saved.values() if set(r["request"]["targets"]) == ready]
        assert len(roots) == 1, "Missing or conflicting focused small-variant acquisition"
        acquired = set()
        def visit(path, receipt):
            pin(path)
            assert receipt["request"]["source_identity"] == identity(source)
            assert receipt["request"]["selection_id"] == selection_id
            assert receipt["request"]["padding"] == selection["variant_padding"]
            verify_input(receipt)
            if receipt["status"] == "partitioned":
                names = []
                for child in receipt["children"]:
                    matches = [(p, r) for p, r in saved.values() if r["request"]["targets"] == child["targets"]]
                    assert len(matches) == 1
                    p, r = matches[0]
                    assert identity(r) == child["receipt_sha256"]
                    names.extend(child["targets"])
                    visit(p, r)
                assert len(names) == len(set(names))
                assert set(names) == set(receipt["request"]["targets"])
            else:
                names = set(receipt["request"]["targets"])
                assert not acquired & names
                acquired.update(names)
                row = dict(source_id=sid, targets=sorted(names), acquisition_status=receipt["status"],
                           input_limitations=receipt.get("receipt", {}).get("limits", []))
                if "bam" in receipt:
                    count = int(subprocess.run(["samtools", "view", "-c", receipt["bam"]["path"]],
                                               capture_output=True, check=True).stdout)
                    row.update(original_records=count, bam_sha256=receipt["bam"]["sha256"],
                               bam_path=str(Path(receipt["bam"]["path"]).relative_to(directory)))
                rows.append(row)
        visit(*roots[0])
        assert acquired == ready
    batches = {gid: b for b in read_json(directory / "batches.json") for gid in b["geometries"]}
    views = []
    engine = None
    for sid in selection["source_ids"]:
        for gid in selection["geometries"]:
            leaves, receipts = acquisition.acquisition_tree(directory, sid, batches[gid]["id"])
            assert len(leaves) == 1
            receipt = next(iter(leaves.values()))
            verify_input(receipt)
            for bid in receipts:
                pin(directory / "sources" / sid / (bid + ".json"))
            native_total = None
            if "bam" in receipt:
                native_total = int(subprocess.run(["samtools", "view", "-c", receipt["bam"]["path"]],
                                                  capture_output=True, check=True).stdout)
            for view in ("forward", "reverse"):
                path = reconstruct.result_path(directory, sid, gid, view)
                result = read_json(path)
                pin(path)
                req = result["request"]
                assert req["acquisition_sha256"] == identity(receipt)
                assert req["inventory_sha256"] == pins["inventory.json.gz"]
                assert req["references_sha256"] == pins["references.json.gz"]
                assert req["source_id"] == sid and req["geometry_id"] == gid
                assert req["orientation"] == view
                assert req["parameters"]["dense_support"]
                if engine is not None:
                    assert engine == req["engine_id"]
                engine = req["engine_id"]
                views.append(dict(source_id=sid, geometry_id=gid, orientation=view, status=result["status"],
                                  candidates=len(result.get("candidates", [])),
                                  limitations=result.get("limitations", []),
                                  input_limitations=result.get("input_limitations", []),
                                  acquisition_status=receipt["status"],
                                  bam_path=(str(Path(receipt["bam"]["path"]).relative_to(directory))
                                            if "bam" in receipt else None),
                                  original_input_records=native_total,
                                  discovery=result.get("discovery"), support_acquisition=result.get("support_acquisition")))
                if "reconstruction" in result:
                    support = result["support_acquisition"]
                    assert support["records_complete"] and support["records_scanned"] == native_total
                    details = directory / result["reconstruction"]["path"]
                    assert digest(details) == result["reconstruction"]["sha256"]
                    pin(details)
                candidates.extend(dict(c, source_id=sid, geometry_id=gid, orientation=view)
                                  for c in result.get("candidates", []))
    for gid in selection["geometries"]:
        available = [c for c in candidates if c["geometry_id"] == gid and c["rna_support"]["fragments"] > 0]
        if not available:
            continue
        c = max(available, key=lambda c: c["rna_support"]["fragments"])
        details = read_json(directory / "reconstructions" / c["source_id"] / (gid + "." + c["orientation"] + ".json.gz"))
        check = check_candidate(c, details, reference)
        receipt = next(iter(acquisition.acquisition_tree(directory, c["source_id"], batches[gid]["id"])[0].values()))
        names = {details["original_records"][rid].split("\t", 1)[0] for rid in check["original_record_hashes"]}
        expected = set(check["original_record_hashes"].values())
        observed = set()
        with pysam.AlignmentFile(receipt["bam"]["path"]) as bam:
            for read in bam:
                if read.query_name in names:
                    observed.add(sha256(read.to_string().encode()).hexdigest())
        assert expected <= observed
        check["original_bam_records_verified"] = len(expected)
        checks.append(check)
    expected_views = 2 * len(selection["geometries"]) * len(selection["source_ids"])
    assert len(views) == expected_views
    output.mkdir(parents=True, exist_ok=True)
    manifest_path = output / "input-files.tsv"
    with manifest_path.open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(["source_id", "source_key", "target_type", "targets", "status", "original_records", "bam_path"])
        for row in rows:
            writer.writerow([row["source_id"], manifest["sources"][row["source_id"]]["key"], "small_variants",
                             ";".join(row["targets"]), row["acquisition_status"],
                             row.get("original_records", ""), row.get("bam_path", "")])
        for row in views:
            if row["orientation"] == "forward":
                writer.writerow([row["source_id"], manifest["sources"][row["source_id"]]["key"], "fusion_geometry",
                                 row["geometry_id"], row["acquisition_status"],
                                 row["original_input_records"], row["bam_path"]])
    ledger = dict(scope=selection["scope"], selection_sha256=digest(directory / "priority-selection.json"),
                  selection=selection, selected_variant_product_pairs=len(ready) * len(selection["source_ids"]),
                  selected_geometry_product_pairs=expected_views // 2, intended_fusion_views=expected_views,
                  all_selected_pairs_accounted=True, full_catalogue_completed=False,
                  small_variant_inputs=rows, fusion_views=views,
                  fusion_outcome_counts=dict(Counter(v["status"] for v in views)),
                  candidates=len(candidates), independent_checks=checks, engine_id=engine, pins=pins,
                  input_files_manifest=dict(path=manifest_path.name, sha256=digest(manifest_path)),
                  pinned_logical_bytes=sum((directory / name).stat().st_size for name in pins))
    write_json(output / "priority-subsets.json.gz", ledger)
    print("Focused subset accounted:", ledger["selected_variant_product_pairs"], "variant/product pairs,",
          expected_views, "fusion views; pinned bytes:", ledger["pinned_logical_bytes"])


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("directory", type=Path)
    parser.add_argument("output", type=Path)
    args = parser.parse_args()
    verify(args.directory, args.output)


if __name__ == "__main__":
    main()
