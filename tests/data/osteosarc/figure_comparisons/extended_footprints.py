"""Selected additional GRCh38 products, never pooled as independent replicates."""

import argparse
from concurrent.futures import ThreadPoolExecutor
import gzip
import json
from pathlib import Path
import shutil

from tests.data.osteosarc.expansion.inventory import digest, fetch_snapshot
from tests.data.osteosarc.figure_comparisons import footprints
from tests.data.osteosarc.figure_comparisons.nr2f2_libraries import load as load_libraries

CORPUS = Path(__file__).parent / "corpus"
PRODUCTS = [
    ("f30f618fb76a0e49", "T0-BostonGene", "T0", "Illumina bulk", "BG003082"),
    ("a5079ad06c1ef4ea", "T0-Personalis", "T0", "Illumina bulk", "Personalis-T0"),
    ("c75479dce2b07bbf", "T1-Tempus", "T1", "Illumina bulk", "Tempus-T1"),
    ("0066232879babe83", "T1-PacBio", "T1", "PacBio", "T1-PacBio"),
    ("62a9d3e67e1b564e", "T3-ONT-dedup", "T3", "ONT", "T3-ONT"),
    ("58c4d68e5f7e55b5", "T3-ONT-tagged", "T3", "ONT", "T3-ONT"),
    ("31e387650e075de2", "T3-scRNA", "T3", "Illumina scRNA", "T3-CellRanger"),
    ("8b8a9a02cbcbcf4b", "T3-CD45neg", "T3", "Illumina scRNA", "T3-CD45neg-CellRanger"),
]


def source_products():
    """Family names track processing provenance, not molecular independence."""
    libraries = {s["source_id"]: s for s in load_libraries()["sources"]}
    sources = []
    for s in footprints.source_products():
        ont = s["product"].startswith("ONT")
        sources.append(dict(s, technology="ONT" if ont else "Illumina bulk",
                            processing_family=(s["sample"] + "-ONT" if ont else
                                               "BG009368" if s["sample"] == "T1" else "SARC0277")))
    for sid, label, sample, technology, family in PRODUCTS:
        original = libraries[sid]
        if original.get("assembly") != "GRCh38" or original["status"] != "ok":
            raise ValueError("Additional source is not a verified GRCh38 product: " + sid)
        product = "ONT-tagged" if label.endswith("tagged") else "ONT-dedup" if label.endswith("dedup") else technology
        sources.append(dict(id=label, source_id=sid, sample=sample, technology=technology,
                            product=product, processing_family=family, url=original["source_url"],
                            sample_metadata=original.get("sample_metadata")))
    return sources


def acquire(output, baseline):
    """Preserve the six original queries; acquire only eight additional products."""
    output.mkdir(parents=True, exist_ok=False)
    old_ids = {s["id"] for s in footprints.source_products()}
    sources = source_products()
    for source in sources:
        if source["id"] not in old_ids:
            continue
        directory = baseline / source["id"]
        receipt = json.loads((directory / "receipt.json").read_text())
        if (receipt["url"] != source["url"] or
                receipt["bam_sha256"] != digest(directory / "regions.bam") or
                receipt["index_sha256"] != digest(directory / "regions.bam.bai")):
            raise ValueError("Original acquisition identity/hash mismatch: " + source["id"])
        shutil.copytree(directory, output / source["id"])

    def fetch(source):
        # Index downloads are pinned explicitly instead of relying on htslib's
        # implicit index cache. Legacy fusion examples are outside this query.
        index = output / (source["id"] + ".source.bai")
        source["index_receipt"] = fetch_snapshot(source["url"] + ".bai", index)
        return footprints.acquire(source, output,
                                  query_regions=footprints.regions()[:-len(footprints.LEGACY)],
                                  index_path=index.resolve())

    with ThreadPoolExecutor(max_workers=3) as pool:
        list(pool.map(fetch, [s for s in sources if s["id"] not in old_ids]))
    (output / "sources.json").write_text(json.dumps(sources, indent=2) + "\n")


def load():
    manifest = json.loads((CORPUS / "extended-footprints-manifest.json").read_text())
    path = CORPUS / manifest["file"]
    if digest(path) != manifest["sha256"]:
        raise ValueError("Extended footprint evidence checksum mismatch")
    return json.loads(gzip.decompress(path.read_bytes()))


def pin(audit_directory, name):
    """Install a completed generated evidence pin, never overwrite a fixture."""
    data = json.loads((audit_directory / "evidence.json").read_text())
    raw = gzip.compress(json.dumps(data, sort_keys=True).encode(), mtime=0)
    path = CORPUS / (name + ".json.gz")
    with path.open("xb") as handle:
        handle.write(raw)
    with (CORPUS / (name + "-manifest.json")).open("x") as handle:
        json.dump(dict(file=path.name, sha256=digest(path)), handle, indent=2)
        handle.write("\n")


def main():
    import logging
    logging.disable(logging.INFO)
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--acquire", type=Path)
    parser.add_argument("--baseline-inputs", type=Path)
    parser.add_argument("--inputs", type=Path)
    parser.add_argument("--output", type=Path)
    parser.add_argument("--pin", type=Path, help="Completed audit directory to pin offline")
    parser.add_argument("--pin-name", choices=("extended-footprints", "dlg5"))
    args = parser.parse_args()
    if args.pin:
        if not args.pin_name:
            parser.error("--pin requires --pin-name")
        pin(args.pin, args.pin_name)
    elif args.acquire:
        if args.baseline_inputs is None:
            parser.error("--acquire requires --baseline-inputs")
        acquire(args.acquire, args.baseline_inputs)
    else:
        if args.inputs is None or args.output is None:
            parser.error("--inputs and --output are required for analysis")
        sources = json.loads((args.inputs / "sources.json").read_text())
        footprints.audit(args.inputs, args.output, sources=sources, all_proteins=True)


if __name__ == "__main__":
    main()
