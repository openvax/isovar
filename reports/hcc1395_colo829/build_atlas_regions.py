"""Build merged GRCh38 intervals for every COLO829 call and SV breakpoint."""

import argparse
from collections import defaultdict
import csv
import gzip
from pathlib import Path
import re


HERE = Path(__file__).parent
SV_FIXTURES = HERE.parents[1] / "tests/data/colo829_whole_run"
BREAKEND = re.compile(r"[\[\]]([^:\[\]]+):(\d+)[\[\]]")


def breakpoints(row):
    """Return one-based breakpoint coordinates from a source VCF record."""
    chrom, pos, _, _, alt, _, _, info = row[:8]
    tags = dict(item.split("=", 1) for item in info.split(";") if "=" in item)
    partner = BREAKEND.search(alt)
    if partner:
        return ((chrom, int(pos)), (partner.group(1), int(partner.group(2))))
    if "END" in tags and tags["END"].isdigit():
        return ((chrom, int(pos)), (tags.get("CHR2", chrom), int(tags["END"])))
    return ((chrom, int(pos)),)


def regions(inventory, sv_vcfs, flank=500):
    """Yield half-open BED regions covering small alleles and SV breakends."""
    with gzip.open(inventory, "rt", newline="") as source:
        for row in csv.DictReader(source, delimiter="\t"):
            if row["sample"] != "COLO829" or row["source"] != "SMaHT_COLO829BLT50":
                continue
            start = max(0, int(row["pos1"]) - 1)
            yield row["chrom1"], start, start + max(1, len(row["ref"]))
    for vcf in sv_vcfs:
        with vcf.open() as source:
            for line in source:
                if line.startswith("#"):
                    continue
                for chrom, pos in breakpoints(line.rstrip("\n").split("\t")):
                    yield chrom if chrom.startswith("chr") else f"chr{chrom}", \
                        max(0, pos - flank - 1), pos + flank


def merged_regions(intervals):
    """Merge overlaps and adjacent intervals, keeping chromosome identity."""
    by_chrom = defaultdict(list)
    for chrom, start, end in intervals:
        by_chrom[chrom].append((start, end))
    for chrom in sorted(by_chrom):
        current_start = current_end = None
        for start, end in sorted(by_chrom[chrom]):
            if current_start is None:
                current_start, current_end = start, end
            elif start <= current_end:
                current_end = max(current_end, end)
            else:
                yield chrom, current_start, current_end
                current_start, current_end = start, end
        if current_start is not None:
            yield chrom, current_start, current_end


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--inventory", type=Path, default=HERE / "somatic_call_inventory.tsv.gz")
    parser.add_argument("--sv-vcf", type=Path, nargs=2,
                        default=(SV_FIXTURES / "colo829-sv-grch38.vcf",
                                 SV_FIXTURES / "colo829-sv-2022-hg38lifted.vcf"))
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    with args.output.open("w") as target:
        for chrom, start, end in merged_regions(regions(args.inventory, args.sv_vcf)):
            target.write(f"{chrom}\t{start}\t{end}\n")


if __name__ == "__main__":
    main()
