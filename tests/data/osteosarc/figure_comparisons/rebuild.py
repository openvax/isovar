"""Pin complete PIP5K1A regions from two public T1 RNA products.

python -m tests.data.osteosarc.figure_comparisons.rebuild ACQUISITION_DIR NEW_OUTPUT_DIR
Uses expansion.acquire receipts; never rewrites an alignment or downsamples.
"""

import argparse
import json
from pathlib import Path

import pysam

from tests.data.osteosarc.expansion.inventory import digest


def rebuild(acquisition, output):
    acquisition, output = Path(acquisition), Path(output)
    output.mkdir(parents=True, exist_ok=False)
    corpus = Path(__file__).resolve().parents[1] / "expansion/corpus"
    variant = next(c["variant"] for c in json.loads((corpus / "manifest.json").read_text())["cases"]
                   if c["variant"]["gene"] == "PIP5K1A")
    cases = []
    for source_id, label in (("53f498a544883d51", "PIP5K1A-T1-ONT"),
                             ("c89442609fffd3f1", "PIP5K1A-T1-Illumina")):
        directory = acquisition / "alignments" / source_id
        receipt = json.loads((directory / "regions-GRCh38.json").read_text())
        source = directory / "regions-GRCh38.bam"
        if receipt["status"] != "ok" or digest(source) != receipt["bam_sha256"]:
            raise ValueError("Acquisition checksum/status mismatch")
        if receipt["variant_ids"] != [variant["variant_id"]]:
            raise ValueError("Expected only the requested PIP5K1A region")
        bam_name = label + ".primary.bam"
        original_hashes, retained = [], []
        from hashlib import sha256
        with pysam.AlignmentFile(source) as bam, pysam.AlignmentFile(output / bam_name, "wb", template=bam) as target:
            for read in bam:
                checksum = sha256(read.to_string().encode()).hexdigest()
                original_hashes.append(checksum)
                if not read.flag & (256 | 1024 | 2048):
                    target.write(read)
                    retained.append(checksum)
        pysam.index(str(output / bam_name))
        cases.append(dict(case_id=label, source_id=source_id, source_url=receipt["request"]["source_url"],
                          variant=variant, primary_bam=bam_name, acquisition=receipt,
                          original_record_sha256=original_hashes, retained_record_sha256=retained,
                          excluded_flags=[256, 1024, 2048],
                          files={name: digest(output / name) for name in (bam_name, bam_name + ".bai")}))
    (output / "manifest.json").write_text(json.dumps(dict(cases=cases), indent=2, sort_keys=True) + "\n")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("acquisition")
    parser.add_argument("output")
    args = parser.parse_args()
    rebuild(args.acquisition, args.output)
