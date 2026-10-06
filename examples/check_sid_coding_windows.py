"""Independently locate selected translated windows in original SAM SEQ/QUAL."""

import argparse
from hashlib import sha256
from pathlib import Path

from Bio.Seq import Seq
import pysam

from examples.sid_sv_audit.inventory import digest, read_json, write_json


def check_window(protein, bam_path):
    expected = {w["fragment_id"]: w["nucleotide_sha256"] for w in protein["witnesses"]}
    found, qualified, records = set(), set(), set()
    aa = protein["amino_acids"] + ("*" if protein["ends_with_stop_codon"] else "")
    with pysam.AlignmentFile(bam_path) as bam:
        for read in bam:
            group = read.get_tag("RG") if read.has_tag("RG") else ""
            fid = sha256(repr((group, read.query_name)).encode()).hexdigest()
            if fid not in expected or read.query_sequence is None:
                continue
            sequence = read.query_sequence
            quality = None if read.query_qualities is None else list(read.query_qualities)
            for reverse in (False, True):
                seq = str(Seq(sequence).reverse_complement()) if reverse else sequence
                scores = quality[::-1] if reverse and quality is not None else quality
                for frame in range(3):
                    coding = seq[frame:frame + ((len(seq) - frame) // 3) * 3]
                    translation = str(Seq(coding).translate())
                    offset = translation.find(aa)
                    while offset >= 0:
                        start = frame + offset * 3
                        nt = seq[start:start + len(aa) * 3]
                        if sha256(nt.encode()).hexdigest() == expected[fid]:
                            found.add(fid)
                            records.add(sha256(read.to_string().encode()).hexdigest())
                            if scores is not None and min(scores[start:start + len(nt)]) >= 20:
                                qualified.add(fid)
                        offset = translation.find(aa, offset + 1)
    if found != set(expected):
        raise ValueError("Claimed full-window witnesses absent from original SAM sequence")
    if len(found) != protein["full_window_fragments"]:
        raise ValueError("Full-window fragment count disagrees")
    if len(qualified) != protein["q20_full_window_fragments"]:
        raise ValueError("Full-window quality count disagrees")
    return dict(full_window_fragments=len(found), native_q20_fragments=len(qualified),
                original_record_hashes=sorted(records), bam_sha256=digest(bam_path),
                amino_acids=protein["amino_acids"], mutation_interval=protein["mutation_interval"],
                ends_with_stop_codon=protein["ends_with_stop_codon"])


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("screen", type=Path)
    parser.add_argument("output", type=Path)
    parser.add_argument("--audit", type=Path, action="append", required=True)
    args = parser.parse_args()
    report = read_json(args.screen)
    checks = []
    for target in sorted({r["target"][0] for r in report["candidates"] if r["kind"] == "indel_frameshift"}):
        available = [r for r in report["candidates"] if r["target"] == [target] and
                     r["kind"] == "indel_frameshift" and r["fragments"] >= 2 and
                     r["annotated_protein_comparison"]["novel_il_collapsed_peptides"]]
        for sid in sorted({r["source_id"] for r in available}):
            rows = [r for r in available if r["source_id"] == sid]
            # Check the longest shifted tail and the strongest complete window.
            chosen = [max(rows, key=lambda r: (len(r["mutant_amino_acids"]), r["q20_full_window_fragments"])),
                      max(rows, key=lambda r: (r["fragments"], len(r["mutant_amino_acids"])))]
            for row in {r["amino_acids"]: r for r in chosen}.values():
                saved = read_json(row["result_path"])
                protein = next(p for p in saved["proteins"] if p["amino_acids"] == row["amino_acids"] and
                               p["mutation_interval"] == row["mutation_interval"])
                # The acquisition receipts identify the original regional BAM.
                bam_path = None
                for base in args.audit:
                    for receipt in (base / "small-variants" / sid).glob("*.json"):
                        acquired = read_json(receipt)
                        if acquired.get("bam", {}).get("sha256") == saved["request"]["bam_sha256"]:
                            bam_path = acquired["bam"]["path"]
                if bam_path is None:
                    raise ValueError("Missing pinned original BAM for coding-window check")
                if digest(bam_path) != saved["request"]["bam_sha256"]:
                    raise ValueError("Original BAM changed")
                checks.append(dict(target_id=target, source_id=sid, result_sha256=digest(row["result_path"]),
                                   **check_window(protein, bam_path)))
    write_json(args.output, dict(candidate_screen_sha256=digest(args.screen), checks=checks,
                                code_sha256=digest(__file__)))
    print("Verified", len(checks), "coding windows against original SAM records")


if __name__ == "__main__":
    main()
