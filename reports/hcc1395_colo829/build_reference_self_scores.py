"""Score top local mutant contexts against pinned Ensembl 87 protein 9-mers.

This is an exact reference-proteome screen, not a patient self or MHC model.
Supply the locally installed Ensembl 87 all-protein FASTA; no download occurs.
"""

import argparse
import csv
import gzip
import hashlib
import io
from pathlib import Path

from isovar.evidence_report import (
    read_tsv, reference_ninemer_matches, self_novelty,
)


HERE = Path(__file__).parent
OUTPUT = HERE / "reference_self_9mer_scores.tsv.gz"
EXPECTED_SHA256 = "d873c199d5dbf20869db3f2cbd6ac726259215d91c2782f9f1ddcad4bb117cf9"
FIELDS = ("sample", "call_id", "novel_9mers", "total_9mers", "reference_fasta_sha256")


def build(fasta, output=OUTPUT):
    digest = hashlib.sha256(Path(fasta).read_bytes()).hexdigest()
    if digest != EXPECTED_SHA256:
        raise ValueError(f"Expected Ensembl 87 protein FASTA {EXPECTED_SHA256}, got {digest}")
    top = {}
    for sample, prefix in (("HCC1395", "hcc1395_direct_rna"),
                           ("COLO829", "colo829")):
        for protein in read_tsv(HERE / f"{prefix}_atlas_proteins.tsv.gz"):
            key = sample, protein["call_id"]
            if key not in top or int(protein["rank"]) < int(top[key]["rank"]):
                top[key] = protein
    matches = reference_ninemer_matches(top.values(), fasta)
    output.parent.mkdir(parents=True, exist_ok=True)
    with output.open("wb") as raw, gzip.GzipFile(filename="", fileobj=raw, mode="wb", mtime=0) as zipped:
        with io.TextIOWrapper(zipped, encoding="utf-8", newline="") as target:
            writer = csv.DictWriter(target, FIELDS, delimiter="\t", lineterminator="\n")
            writer.writeheader()
            for (sample, call_id), protein in sorted(top.items()):
                novel, total = self_novelty(protein, matches)
                writer.writerow(dict(sample=sample, call_id=call_id, novel_9mers=novel,
                                     total_9mers=total, reference_fasta_sha256=digest))
    return len(top)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--fasta", type=Path, required=True)
    parser.add_argument("--output", type=Path, default=OUTPUT)
    args = parser.parse_args()
    print("Scored local products:", build(args.fasta, args.output))


if __name__ == "__main__":
    main()
