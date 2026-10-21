"""Recover the pinned COLO829 original SRA pairs used by the RNA fixture."""

import argparse
from collections import defaultdict
import csv
from pathlib import Path
import subprocess


HERE = Path(__file__).parent


def extract(fastq_dump, spot_table, output_dir):
    """Write paired FASTQs for the exact source spot IDs in the audit table."""
    output_dir.mkdir(parents=True, exist_ok=True)
    by_gene = defaultdict(set)
    with spot_table.open(newline="") as source:
        for row in csv.DictReader(source, delimiter="\t"):
            by_gene[row["gene"]].add(int(row["sra_spot_id"]))
    written_ids = set()
    with (output_dir / "R1.fastq").open("w") as r1, \
            (output_dir / "R2.fastq").open("w") as r2:
        for gene, ids in by_gene.items():
            command = [str(fastq_dump), "-N", str(min(ids)), "-X", str(max(ids)),
                       "--split-spot", "-Z", "SRR8615617"]
            with subprocess.Popen(command, stdout=subprocess.PIPE, text=True) as process:
                while True:
                    first = [process.stdout.readline() for _ in range(4)]
                    if not first[0]:
                        break
                    second = [process.stdout.readline() for _ in range(4)]
                    if not second[0]:
                        raise ValueError(f"Incomplete SRA spot near {gene}")
                    spot = int(first[0].split()[0].split(".")[-1])
                    if spot in ids:
                        if spot in written_ids:
                            raise ValueError(f"Duplicate SRA spot {spot}")
                        r1.writelines(first)
                        r2.writelines(second)
                        written_ids.add(spot)
                if process.wait():
                    raise RuntimeError(f"fastq-dump failed at {gene}")
    expected = set().union(*by_gene.values())
    if written_ids != expected:
        raise ValueError(f"Missing SRA spot IDs: {expected - written_ids}")
    return len(written_ids)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--fastq-dump", type=Path, required=True)
    parser.add_argument("--spot-table", type=Path,
                        default=HERE / "colo829_original_spot_ids.tsv")
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()
    print(f"Recovered {extract(args.fastq_dump, args.spot_table, args.output_dir)} SRA pairs")


if __name__ == "__main__":
    main()
