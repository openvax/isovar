"""Manually rebuild the CC0 osteosarc.com read fixtures; never downloads in CI.

Usage: python tests/data/osteosarc/rebuild.py SOURCE_DIR OUTPUT_DIR
SOURCE_DIR contains bulk_star_t0.sam, ont_t1.sam, selection.json, vafs.tsv,
vafs-columns.tsv and data-license.yaml (see README for acquisition).
"""

import argparse
from collections import Counter
import csv
import gzip
from hashlib import sha256
import json
from pathlib import Path
import sys

import pysam

sys.path.insert(0, str(Path(__file__).resolve().parents[3]))
from tests.real_rna_helpers import cigar_observation, record_digest  # noqa: E402


def build(source_dir, output_dir):
    output_dir.mkdir(parents=True, exist_ok=True)
    selection = json.loads((source_dir / "selection.json").read_text())
    variants = [dict(v, pos=int(v["pos"])) for v in selection["variants"]]
    manifest = {
        "schema_version": 1,
        "assembly": "GRCh38",
        "license": "CC0-1.0",
        "accessed": "2026-09-07",
        "selection": (
            "First 48 template names by SHA256(name) per locus, allele/quality blind, "
            "plus two lowest-quality primary alternate-supporting template names per locus; "
            "all available regional records for those names. Not suitable for estimating VAF."
        ),
        "datasets": {},
    }
    for label, source in selection["sources"].items():
        with pysam.AlignmentFile(source_dir / (label + ".sam")) as bam:
            reads = list(bam)
            header = bam.text
        selected_names = set()
        for variant in variants:
            names = {r.query_name for r in reads if r.reference_name == variant["chrom"]
                     and r.get_overlap(variant["pos"] - 1, variant["pos"] + len(variant["ref"])) > 0}
            selected_names.update(sorted(names, key=lambda n: sha256(n.encode("ascii")).hexdigest())[:48])
        baseline_names = sorted(selected_names)
        quality_stress_names = set()
        for variant in variants:
            alt = variant["alt"] if len(variant["ref"]) == 1 else variant["alt"][1:]
            qualities_by_name = {}
            for read in reads:
                if read.flag & (256 | 1024 | 2048):
                    continue
                observation = cigar_observation(read, variant)
                if observation is None or observation["allele"] != alt or not observation["qualities"]:
                    continue
                quality = min(observation["qualities"])
                qualities_by_name[read.query_name] = min(quality, qualities_by_name.get(read.query_name, quality))
            quality_stress_names.update(sorted(qualities_by_name, key=lambda n: (
                qualities_by_name[n], sha256(n.encode("ascii")).hexdigest()))[:2])
        selected_names.update(quality_stress_names)
        selected = [r for r in reads if r.query_name in selected_names]
        raw = (header + "".join(r.to_string() + "\n" for r in selected)).encode("ascii")
        filename = label + ".sam.gz"
        with (output_dir / filename).open("wb") as handle:
            with gzip.GzipFile(filename="", mode="wb", fileobj=handle, mtime=0) as compressed:
                compressed.write(raw)
        cases = []
        for variant in variants:
            observations = [o for r in selected if (o := cigar_observation(r, variant)) is not None]
            full = [o for r in reads if (o := cigar_observation(r, variant)) is not None
                    and not r.flag & (256 | 1024 | 2048)]
            alt = variant["alt"] if len(variant["ref"]) == 1 else variant["alt"][1:]
            # All selected indels here are anchored deletions.
            cases.append(dict(variant, observations=observations,
                              source_primary_observations=len(full),
                              source_primary_alt_quality_retention={str(q): sum(
                                  o["allele"] == alt and o["qualities"] is not None
                                  and min(o["qualities"]) >= q for o in full)
                                  for q in (0, 10, 20, 30)}))
        manifest["datasets"][label] = {
            "source": source,
            "source_sam_sha256": sha256((source_dir / (label + ".sam")).read_bytes()).hexdigest(),
            "source_region_records": len(reads),
            "file": filename,
            "sam_sha256": sha256(raw).hexdigest(),
            "fixture_records": len(selected),
            "baseline_template_names": baseline_names,
            "quality_stress_template_names": sorted(quality_stress_names),
            "records": [record_digest(r) for r in selected],
            "variants": cases,
        }
        print(label, len(selected), [(v["gene"], dict(Counter(o["allele"] for o in v["observations"]))) for v in cases])
    (output_dir / "manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
    # Preserve original source lines for these variants; counts remain source
    # claims, not expectations derived from our downsampled read collection.
    with (source_dir / "vafs.tsv").open() as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        variant_ids = {v["variant_id"] for v in variants}
        rows = [r for r in reader if r["variant_id"] in variant_ids]
        with (output_dir / "source_variant_counts.tsv").open("w") as output:
            writer = csv.DictWriter(output, fieldnames=reader.fieldnames, delimiter="\t", lineterminator="\n")
            writer.writeheader()
            writer.writerows(rows)
    for source, destination in (("vafs-columns.tsv", "source_variant_counts.columns.tsv"),
                                ("data-license.yaml", "source_registry.yaml")):
        (output_dir / destination).write_bytes((source_dir / source).read_bytes())
    (output_dir / "source_checksums.json").write_text(json.dumps({
        name: sha256((source_dir / name).read_bytes()).hexdigest()
        for name in ("vafs.tsv", "vafs-columns.tsv", "data-license.yaml", "bams.json")}, indent=2) + "\n")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("source_dir", type=Path)
    parser.add_argument("output_dir", type=Path)
    args = parser.parse_args()
    build(args.source_dir, args.output_dir)
