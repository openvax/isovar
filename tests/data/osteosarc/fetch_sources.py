"""Optional network acquisition; fetch indexed loci, never entire BAMs or in CI.

Usage: python tests/data/osteosarc/fetch_sources.py NEW_SOURCE_DIR
Requires curl and samtools with HTTPS support. Metadata drift is an error;
revising the pinned source snapshots requires explicit review.
"""

import argparse
from concurrent.futures import ThreadPoolExecutor
from hashlib import sha256
import json
from pathlib import Path
import subprocess


HERE = Path(__file__).resolve().parent
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
        subprocess.run(["curl", "--fail", "--location", "--retry", "3", "--output", str(path), url], check=True)
        if sha256(path.read_bytes()).hexdigest() != checksums[filename]:
            raise ValueError(f"Upstream metadata changed: {filename}; review before rebuilding")
    regions = [f"{v['chrom']}:{int(v['pos']) - 1}-{int(v['pos']) + len(v['ref'])}"
               for v in selection["variants"]]

    def fetch_bam(item):
        label, source = item
        # -M avoids duplicate emission when intervals overlap. Original
        # duplicate SAM records already in a source BAM remain untouched.
        subprocess.run(["samtools", "view", "--no-PG", "-h", "-M", "-o",
                        str(destination / (label + ".sam")), source["url"], *regions],
                       check=True, cwd=destination)
        print(label, "downloaded", flush=True)

    with ThreadPoolExecutor(max_workers=2) as executor:
        list(executor.map(fetch_bam, selection["sources"].items()))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("destination", type=Path)
    fetch(parser.parse_args().destination)
