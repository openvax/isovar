"""Extract the selected source VCF records for the original RNA fixtures."""

import argparse
import csv
import gzip
from pathlib import Path


HERE = Path(__file__).parent
DEFAULT_OUTPUT = Path(__file__).parents[2] / "tests/data/hcc1395_colo829_original"


def extract(source_path, target_path, selected):
    """Preserve source headers and exact records for selected alleles."""
    found = set()
    with gzip.open(source_path, "rt") as source, target_path.open("w") as target:
        for line in source:
            if line.startswith("#"):
                target.write(line)
                continue
            chrom, pos, _, ref, alt, *_ = line.split("\t")
            key = chrom, pos, ref, alt
            if key in selected:
                if key in found:
                    raise ValueError(f"Duplicate source VCF record: {key}")
                target.write(line)
                found.add(key)
    if found != selected:
        raise ValueError(f"Missing source VCF alleles: {selected - found}")
    return len(found)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--hcc-vcf", type=Path, required=True)
    parser.add_argument("--colo-vcf", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT)
    args = parser.parse_args()
    with (HERE / "original_rna_candidates.tsv").open(newline="") as source:
        candidates = list(csv.DictReader(source, delimiter="\t"))
    args.output_dir.mkdir(parents=True, exist_ok=True)
    for sample, source, name in (
            ("HCC1395", args.hcc_vcf, "hcc1395_selected.vcf"),
            ("COLO829", args.colo_vcf, "colo829_selected.vcf")):
        selected = {(row["chrom"], row["pos"], row["ref"], row["alt"])
                    for row in candidates if row["sample"] == sample}
        count = extract(source, args.output_dir / name, selected)
        print(f"{sample}: {count} source VCF records")


if __name__ == "__main__":
    main()
