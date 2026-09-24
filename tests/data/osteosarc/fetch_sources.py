"""Optional network acquisition; fetch indexed loci, never entire BAMs or in CI.

Usage: python tests/data/osteosarc/fetch_sources.py NEW_SOURCE_DIR
Requires samtools and an explicit ISOVAR_SID_SNAPSHOT. Metadata drift is an error;
revising the pinned source snapshots requires explicit review.
"""

import argparse
from concurrent.futures import ThreadPoolExecutor
from hashlib import sha256
import json
from pathlib import Path
import shutil
import sys

import pysam


HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parents[2]))
from isovar.sid_data import extract_regions, fetch_metadata  # noqa: E402
METADATA = {
    "bams.json": "https://osteosarc.com/bams/bams.json",
    "vafs.tsv": "https://osteosarc.com/variants/variant_vafs_long.tsv",
    "vafs-columns.tsv": "https://osteosarc.com/variants/variant_vafs_long.columns.tsv",
    "data-license.yaml": "https://raw.githubusercontent.com/awslabs/open-data-registry/main/datasets/sid-osteosarc.yaml",
}


def fetch(destination):
    destination = destination.resolve()
    destination.mkdir(parents=True, exist_ok=False)
    checksums = json.loads((HERE / "source_checksums.json").read_text())
    selection = json.loads((HERE / "selection.json").read_text())
    (destination / "selection.json").write_bytes((HERE / "selection.json").read_bytes())
    for filename, url in METADATA.items():
        path = destination / filename
        cached, receipt = fetch_metadata(url)
        shutil.copyfile(cached, path)
        path.with_name(path.name + ".receipt.json").write_text(json.dumps(receipt, indent=2) + "\n")
        if sha256(path.read_bytes()).hexdigest() != checksums[filename]:
            raise ValueError(f"Upstream metadata changed: {filename}; review before rebuilding")
    regions = [f"{v['chrom']}:{int(v['pos']) - 1}-{int(v['pos']) + len(v['ref'])}"
               for v in selection["variants"]]

    def fetch_bam(item):
        label, source = item
        subset = extract_regions(source["url"], regions, "GRCh38")
        with subset.open() as bam, pysam.AlignmentFile(destination / (label + ".sam"), "w", header=bam.header) as sam:
            for read in bam:
                sam.write(read)
        (destination / (label + ".receipt.json")).write_text(json.dumps(subset.receipt, indent=2) + "\n")
        print(label, "downloaded", flush=True)

    with ThreadPoolExecutor(max_workers=2) as executor:
        list(executor.map(fetch_bam, selection["sources"].items()))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("destination", type=Path)
    fetch(parser.parse_args().destination)
