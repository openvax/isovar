"""Rebuild the HCC1395 whole-run original-RNA small-variant atlas."""

import argparse
from pathlib import Path

from .build_colo829_atlas import write_atlas


HERE = Path(__file__).parent


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--inventory", type=Path, default=HERE / "somatic_call_inventory.tsv.gz")
    parser.add_argument("--bam", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, default=HERE)
    parser.add_argument("--rna-run", choices=("SRR23456709", "SRR23456710"),
                        default="SRR23456709")
    parser.add_argument("--prefix")
    args = parser.parse_args()
    prefix = args.prefix or ("hcc1395_direct_rna" if args.rna_run == "SRR23456709"
                             else "hcc1395_cdna")
    print("calls, proteins, ALT witnesses:",
          write_atlas(args.inventory, args.bam, args.output_dir,
                      sample="HCC1395", dna_source="SEQC2_v1.2.1",
                      rna_run=args.rna_run, prefix=prefix))


if __name__ == "__main__":
    main()
