"""Account for every nominated target/product, including incomplete work."""

from collections import Counter
import csv
from pathlib import Path

from .acquisition import acquisition_tree, selected_sources
from .inventory import digest, identity, load_inventory, read_json, write_json
from .reconstruct import result_path


def coverage(manifest, sources, outcomes, headers):
    """No missing job can become a negative RNA call or a completed audit."""
    expected = {(sid, gid, orientation) for sid in sources for gid in manifest["geometries"]
                for orientation in ("forward", "reverse")}
    unexpected = set(outcomes) - expected
    if unexpected:
        raise ValueError("Results outside the intended event/source matrix")
    rows = []
    for sid in sorted(sources):
        for name, target in sorted(manifest["nominations"].items()):
            gid = target["geometry_id"]
            row = dict(source_id=sid, nomination=name, geometry_id=gid)
            if gid is None:
                row.update(status="unresolved_geometry", views={})
            elif sid not in headers:
                row.update(status="pending_header", views={})
            elif headers[sid]["status"] != "ready":
                row.update(status="not_assessable", reason=headers[sid]["status"], views={})
            else:
                views = {view: outcomes.get((sid, gid, view), {}).get("status", "pending_reconstruction")
                         for view in ("forward", "reverse")}
                row.update(status="pending_reconstruction" if "pending_reconstruction" in views.values()
                           else "accounted", views=views)
            rows.append(row)
    counts = Counter(row["status"] for row in rows)
    return dict(expected_nominations=len(manifest["nominations"]), expected_sources=len(sources),
                expected_pairs=len(manifest["nominations"]) * len(sources),
                expected_geometry_views=len(expected), observed_geometry_views=len(outcomes),
                all_pairs_accounted=not any(key.startswith("pending") for key in counts),
                counts=dict(counts), rows=rows)


def report(directory, output, cohort="tumor_candidate", require_complete=False):
    directory, output = Path(directory), Path(output)
    manifest = load_inventory(directory)
    inventory_sha256 = digest(directory / "inventory.json.gz")
    references_sha256 = digest(directory / "references.json.gz")
    sources = selected_sources(manifest, cohort)
    outcomes, headers, candidates, pins = {}, {}, [], {}
    run_identity = None
    batches = read_json(directory / "batches.json")
    batch_definitions = {batch["id"]: batch for batch in batches}
    batch_for_geometry = {}
    for batch in batches:
        for gid in batch["geometries"]:
            if gid in batch_for_geometry:
                raise ValueError("Geometry occurs in more than one acquisition batch")
            batch_for_geometry[gid] = batch["id"]
    if set(batch_for_geometry) != set(manifest["geometries"]):
        raise ValueError("Acquisition batches do not cover the intended geometries")
    expected_filenames = {gid + "." + view + ".json.gz" for gid in manifest["geometries"]
                          for view in ("forward", "reverse")}
    for sid in sorted(sources):
        actual_filenames = {p.name for p in (directory / "results" / sid).glob("*.json.gz")}
        if actual_filenames - expected_filenames:
            raise ValueError("Unexpected result files for selected source")
        acquisition_pins = {}
        acquisition_for_geometry = {}
        hp = directory / "sources" / sid / "header.json"
        if hp.exists():
            headers[sid] = read_json(hp)
            if headers[sid]["source_id"] != sid or headers[sid]["source_identity"] != identity(sources[sid]):
                raise ValueError("Header belongs to another source")
            pins[str(hp.relative_to(directory))] = digest(hp)
        for gid in manifest["geometries"]:
            for view in ("forward", "reverse"):
                path = result_path(directory, sid, gid, view)
                if not path.exists():
                    continue
                result = read_json(path)
                req = result["request"]
                if (req["source_id"], req["geometry_id"], req["orientation"]) != (sid, gid, view):
                    raise ValueError("Result is filed under the wrong identity")
                if req["inventory_sha256"] != inventory_sha256:
                    raise ValueError("Result belongs to another inventory")
                if req["references_sha256"] != references_sha256:
                    raise ValueError("Result belongs to another reference snapshot")
                bid = batch_for_geometry[gid]
                if bid not in acquisition_pins:
                    leaves, receipts = acquisition_tree(directory, sid, bid)
                    if receipts[bid]["request"]["batch"] != batch_definitions[bid]:
                        raise ValueError("Acquisition batch differs from the pinned plan")
                    for aid, acquired in receipts.items():
                        if acquired["request"]["source_identity"] != identity(sources[sid]):
                            raise ValueError("Acquisition belongs to another source")
                        ap = directory / "sources" / sid / (aid + ".json")
                        pins[str(ap.relative_to(directory))] = digest(ap)
                    acquisition_pins[bid] = identity(receipts[bid])
                    for acquired in leaves.values():
                        leaf_sha256 = identity(acquired)
                        for name in acquired["request"]["batch"]["geometries"]:
                            acquisition_for_geometry[name] = leaf_sha256
                if req["acquisition_sha256"] != acquisition_for_geometry[gid]:
                    raise ValueError("Acquisition changed after reconstruction")
                run = {key: req[key] for key in ("parameters", "isovar_version", "engine_id", "max_seconds")}
                if run_identity is not None and run != run_identity:
                    raise ValueError("Results mix engine versions or reconstruction parameters")
                run_identity = run
                if "reconstruction" in result:
                    details = directory / "reconstructions" / sid / path.name
                    if result["reconstruction"]["path"] != str(details.relative_to(directory)):
                        raise ValueError("Reconstruction is filed under the wrong identity")
                    if digest(details) != result["reconstruction"]["sha256"]:
                        raise ValueError("Detailed reconstruction checksum mismatch")
                    pins[str(details.relative_to(directory))] = digest(details)
                outcomes[sid, gid, view] = dict(status=result["status"],
                    limitations=result.get("limitations", []), acquisition_status=result["acquisition_status"],
                    input_limitations=result.get("input_limitations", []), reason=result.get("reason"),
                    orf_search_truncated=result.get("orf_search_truncated"),
                    reference_exclusions=result.get("reference_exclusions", []),
                    unsupported_reference_cds=result.get("unsupported_reference_cds", []))
                pins[str(path.relative_to(directory))] = digest(path)
                for candidate in result.get("candidates", []):
                    candidates.append(dict(source_id=sid, source_url=sources[sid]["url"], geometry_id=gid,
                                           orientation=view, acquisition_status=result["acquisition_status"],
                                           input_limitations=result.get("input_limitations", []), **candidate))
    ledger = coverage(manifest, sources, outcomes, headers)
    ledger.update(inventory_sha256=inventory_sha256, references_sha256=references_sha256,
                  run_identity=run_identity, cohort=cohort,
                  outcome_counts=dict(Counter(row["status"] for row in outcomes.values())),
                  outcomes=[dict(source_id=s, geometry_id=g, orientation=v, **row)
                            for (s, g, v), row in sorted(outcomes.items())])
    write_json(output / "coverage.json.gz", ledger)
    write_json(output / "candidate-orfs.json.gz", candidates)
    write_json(output / "result-pins.json.gz", pins)
    csv_path = output / "candidate-orfs.csv"
    if csv_path.exists():
        raise FileExistsError(csv_path)
    fields = ("source_id", "source_url", "geometry_id", "orientation", "amino_acids", "ends_with_stop_codon",
              "full_interval_fragments", "native_query_fragments", "reverse_complement_fragments",
              "missing_quality_reads", "acquisition_status", "input_limitations", "uncertainty_flags")
    with csv_path.open("x", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        for c in sorted(candidates, key=lambda c: (-c["rna_support"]["fragments"], c["source_id"], c["candidate_id"])):
            support = c["rna_support"]
            row = {key: c[key] for key in fields[:6]}
            row.update(full_interval_fragments=support["fragments"],
                       native_query_fragments=support["fragment_query_orientations"]["original_query"],
                       reverse_complement_fragments=support["fragment_query_orientations"]["reverse_complement"],
                       missing_quality_reads=support["missing_quality_reads"],
                       acquisition_status=c["acquisition_status"],
                       input_limitations=";".join(c["input_limitations"]),
                       uncertainty_flags=";".join(c["uncertainty_flags"]))
            writer.writerow(row)
    if require_complete and not ledger["all_pairs_accounted"]:
        raise ValueError("Audit is incomplete: see pending entries in coverage.json.gz")
    return {key: ledger[key] for key in ("expected_pairs", "all_pairs_accounted", "counts", "outcome_counts")}
