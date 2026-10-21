"""Build GRCh38 intervals for all SEQC2 calls and Talsania SV ends."""

import argparse
import csv
import gzip
from pathlib import Path

from .build_atlas_regions import merged_regions


HERE = Path(__file__).parent


def regions(inventory, flank=500):
    """Yield half-open BED intervals for the HCC1395 source inventory."""
    with gzip.open(inventory, "rt", newline="") as source:
        for row in csv.DictReader(source, delimiter="\t"):
            if row["sample"] != "HCC1395":
                continue
            if row["source"] == "SEQC2_v1.2.1" and row["class"] in ("SNV", "indel"):
                start = max(0, int(row["pos1"]) - 1)
                yield row["chrom1"], start, start + max(1, len(row["ref"]))
            elif row["source"] == "Talsania_2022_S4":
                for side in ("1", "2"):
                    chrom, pos = row[f"chrom{side}"], row[f"pos{side}"]
                    if chrom and pos:
                        pos = int(pos)
                        yield chrom, max(0, pos - flank - 1), pos + flank


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--inventory", type=Path, default=HERE / "somatic_call_inventory.tsv.gz")
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    with args.output.open("w") as target:
        for chrom, start, end in merged_regions(regions(args.inventory)):
            target.write(f"{chrom}\t{start}\t{end}\n")


if __name__ == "__main__":
    main()
