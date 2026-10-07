"""Verify the complete selected recovery scope against original BAM records.

Native samtools recounts every input and its requested regional records.
Biopython translates every ORF; the best-supported candidate per pair also
gets an independent SAM sequence/fragment recount. This does not establish
translation, tumor specificity, RNA polarity or independent cell abundance.
"""

import argparse
from collections import Counter
from hashlib import sha256
from pathlib import Path
import subprocess
import tempfile

from Bio.Seq import Seq
import pysam

from examples.check_sid_sv_pilot import check_candidate
from examples.recover_sid_sv_pilot import coverage
from examples.sid_sv_audit import acquisition, reconstruct, references
from examples.sid_sv_audit.inventory import digest, identity, load_inventory, read_json, write_json


def verify(directory, output):
    """Require all selected reconstructions and verify counts without pooling."""
    directory, output = Path(directory), Path(output)
    manifest, reference = load_inventory(directory), references.load(directory)
    selection = read_json(directory / "recovery-selection.json")
    run = read_json(directory / "recovery-run.json")
    batches = {b["geometries"][0]: b for b in read_json(directory / "batches.json")}
    observed = coverage(directory, selection, batches)
    if not observed["all_inputs_acquired"] or not observed["all_views_reconstructed"]:
        raise ValueError("Recovery verification requires every selected input and reconstruction")
    pins = {name: digest(directory / name) for name in (
        "inventory.json.gz", "references.json.gz", "recovery-selection.json", "recovery-run.json", "batches.json")}
    rows, checks, outcomes = [], [], Counter()
    translations = 0
    for pair in selection["pairs"]:
        sid, gid = pair["source_id"], pair["geometry_id"]
        batch = batches[gid]
        leaves, receipts = acquisition.acquisition_tree(directory, sid, batch["id"])
        assert len(leaves) == 1
        receipt = next(iter(leaves.values()))
        assert receipt["status"] in acquisition.USABLE_STATUSES
        assert receipt["request"]["source_identity"] == identity(manifest["sources"][sid])
        for key in ("bam", "index"):
            assert digest(receipt[key]["path"]) == receipt[key]["sha256"]
        for bid in receipts:
            path = directory / "sources" / sid / (bid + ".json")
            pins[str(path.relative_to(directory))] = digest(path)
        bam_path = receipt["bam"]["path"]
        native_total = int(subprocess.run(["samtools", "view", "-c", bam_path],
                                         capture_output=True, check=True).stdout)
        assert native_total == receipt["receipt"]["records"]
        with tempfile.TemporaryDirectory() as temporary:
            bed = Path(temporary) / "regions.bed"
            bed.write_text("".join("%s\t%d\t%d\n" % tuple(region) for region in batch["regions"]))
            native_regional = int(subprocess.run(["samtools", "view", "-c", "-M", "-L", str(bed), bam_path],
                                                capture_output=True, check=True).stdout)
        candidates, views = [], {}
        for orientation in ("forward", "reverse"):
            path = reconstruct.result_path(directory, sid, gid, orientation)
            result = read_json(path)
            request = result["request"]
            assert request["engine_id"] == run["engine_id"]
            assert request["parameters"]["dense_support"]
            assert request["acquisition_sha256"] == identity(receipt)
            assert request["inventory_sha256"] == pins["inventory.json.gz"]
            assert request["references_sha256"] == pins["references.json.gz"]
            assert (request["source_id"], request["geometry_id"], request["orientation"]) == (sid, gid, orientation)
            support = result["support_acquisition"]
            assert support["records_complete"] and support["records_scanned"] == native_total
            detail_path = directory / result["reconstruction"]["path"]
            assert digest(detail_path) == result["reconstruction"]["sha256"]
            pins.update({str(p.relative_to(directory)): digest(p) for p in (path, detail_path)})
            views[orientation] = {key: result[key] for key in (
                "status", "records", "discovery", "support_acquisition", "input_limitations")}
            views[orientation]["candidates"] = len(result["candidates"])
            outcomes[result["status"]] += 1
            for candidate in result["candidates"]:
                nt = candidate["nucleotide_sequence"]
                aa = str(Seq(nt[:len(nt) // 3 * 3]).translate(table=1))
                assert aa.rstrip("*") == candidate["amino_acids"]
                assert aa.endswith("*") == candidate["ends_with_stop_codon"]
                translations += 1
                candidates.append(dict(candidate, source_id=sid, geometry_id=gid, orientation=orientation))
        if candidates:
            candidate = max(candidates, key=lambda c: c["rna_support"]["fragments"])
            result = read_json(reconstruct.result_path(directory, sid, gid, candidate["orientation"]))
            details = read_json(directory / result["reconstruction"]["path"])
            check = check_candidate(candidate, details, reference)
            expected = set(check["original_record_hashes"].values())
            names = {details["original_records"][key].split("\t", 1)[0]
                     for key in check["original_record_hashes"]}
            found = set()
            with pysam.AlignmentFile(bam_path) as bam:
                for read in bam:
                    if read.query_name in names:
                        found.add(sha256(read.to_string().encode()).hexdigest())
            assert expected <= found
            check["original_bam_records_verified"] = len(expected)
            checks.append(check)
        rows.append(dict(pair=pair["pair"], source_id=sid, geometry_id=gid,
                         nominations=manifest["geometries"][gid]["nominations"],
                         acquisition_status=receipt["status"], regional_alignment_records=native_regional,
                         acquired_records=native_total, bam_sha256=receipt["bam"]["sha256"], views=views))
    ledger = dict(scope=selection["scope"], intended_pairs=observed["expected_pairs"],
                  intended_views=observed["expected_views"], accounted_views=sum(outcomes.values()),
                  outcome_counts=dict(outcomes), rows=rows, independent_checks=checks,
                  independent_orf_translations=translations, pins=pins,
                  code_sha256=digest(__file__), full_catalogue_complete=False)
    assert ledger["accounted_views"] == ledger["intended_views"]
    write_json(output / "recovery-checks.json.gz", ledger)
    return ledger


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("directory", type=Path)
    parser.add_argument("output", type=Path)
    args = parser.parse_args()
    ledger = verify(args.directory, args.output)
    print(ledger["outcome_counts"], "independently translated ORFs:", ledger["independent_orf_translations"])


if __name__ == "__main__":
    main()
