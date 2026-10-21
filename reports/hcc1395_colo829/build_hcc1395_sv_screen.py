"""Screen original HCC1395 ONT RNA near every Talsania SV breakend."""

import argparse
import csv
import gzip
from pathlib import Path

import pysam

from .build_colo829_sv_screen import split_candidates, window_alignments


HERE = Path(__file__).parent
FIELDS = (
    "sample", "dna_source", "class", "call_id", "chrom1", "pos1", "chrom2",
    "pos2", "dna_filter", "dna_support", "published_rna_fusion",
    "rna_run", "window_bp", "end1_alignment_records", "end2_alignment_records",
    "end1_reads", "end2_reads", "observed_split_alignment_candidates",
    "candidate_read_ids", "rna_status", "path_status", "protein_status",
)


def source_sv_rows(inventory):
    """Yield Talsania SV records without collapsing DNA event identities."""
    with gzip.open(inventory, "rt", newline="") as source:
        for row in csv.DictReader(source, delimiter="\t"):
            if row["sample"] == "HCC1395" and row["source"] == "Talsania_2022_S4":
                yield row


def screen_rows(inventory, alignment, flank=500, rna_run="SRR23456709"):
    """Yield an RNA coverage and split-alignment screen for each SV record."""
    for source in source_sv_rows(inventory):
        chrom1, pos1 = source["chrom1"], int(source["pos1"])
        chrom2, pos2 = source["chrom2"], int(source["pos2"])
        first = window_alignments(alignment, chrom1, pos1, flank)
        second = window_alignments(alignment, chrom2, pos2, flank)
        distinct_windows = chrom1 != chrom2 or abs(pos1 - pos2) > 2 * flank
        witnesses = split_candidates(first, second) if distinct_windows else []
        if witnesses:
            status = "split_alignment_candidate_unresolved"
        elif first and second:
            status = "both_ends_covered_no_split_candidate"
        elif first or second:
            status = "one_end_covered_no_split_candidate"
        else:
            status = "no_breakend_window_alignment"
        yield {
            "sample": "HCC1395", "dna_source": source["source"],
            "class": source["class"], "call_id": source["call_id"],
            "chrom1": chrom1, "pos1": pos1, "chrom2": chrom2, "pos2": pos2,
            "dna_filter": source["filter"], "dna_support": source["dna_support"],
            "published_rna_fusion": source["published_rna_fusion"],
            "rna_run": rna_run, "window_bp": flank,
            "end1_alignment_records": len(first),
            "end2_alignment_records": len(second),
            "end1_reads": len({r.query_name for r in first}),
            "end2_reads": len({r.query_name for r in second}),
            "observed_split_alignment_candidates": len(witnesses),
            "candidate_read_ids": ";".join(witnesses),
            "rna_status": status, "path_status": "unresolved",
            "protein_status": "not_called_unresolved_path",
        }


def write_screen(inventory, bam, output, flank=500, rna_run="SRR23456709"):
    """Write all SV screen rows; return the number written."""
    with pysam.AlignmentFile(str(bam), "rb") as alignment, output.open("w", newline="") as target:
        writer = csv.DictWriter(target, FIELDS, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        count = 0
        for row in screen_rows(inventory, alignment, flank, rna_run):
            writer.writerow(row)
            count += 1
    return count


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--inventory", type=Path, default=HERE / "somatic_call_inventory.tsv.gz")
    parser.add_argument("--bam", type=Path, required=True)
    parser.add_argument("--output", type=Path)
    parser.add_argument("--rna-run", choices=("SRR23456709", "SRR23456710"),
                        default="SRR23456709")
    args = parser.parse_args()
    default_name = ("hcc1395_direct_rna_sv_screen.tsv" if args.rna_run == "SRR23456709"
                    else "hcc1395_cdna_sv_screen.tsv")
    output = args.output or HERE / default_name
    print(f"Screened {write_screen(args.inventory, args.bam, output, rna_run=args.rna_run)} SV records")


if __name__ == "__main__":
    main()
