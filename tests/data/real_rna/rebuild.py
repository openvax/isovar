"""Rebuild public RNA fixtures from downloaded source regions; never run in CI.

See README.md for download commands and upstream provenance. This script does
not import isovar and cannot silently update expectations from its output.
"""

import argparse
from collections import Counter
import gzip
from hashlib import sha256
import json
from pathlib import Path
import sys

import pysam

sys.path.insert(0, str(Path(__file__).resolve().parents[3]))
from tests.real_rna_helpers import cigar_observation, record_digest  # noqa: E402


CASES = {
    "illumina_star": [("chr4", p) for p in (337869, 337877, 337904, 768474)],
    "ont_minimap2": [("chr6", p) for p in (43770613, 43780233, 43785475, 43785588)],
    "pacbio_minimap2": [],
}
SOURCES = {
    "illumina_star": "star.bam",
    "ont_minimap2": "ont_chr6.sam",
    "pacbio_minimap2": "pacbio_k562.sam",
}


def build(source_dir, output_dir):
    output_dir.mkdir(parents=True, exist_ok=True)
    variants, vcf_headers = {}, []
    for chromosome in ("chr4", "chr6"):
        for line in (source_dir / (chromosome + ".vcf")).read_text().splitlines():
            if line.startswith("#"):
                if chromosome == "chr4":
                    vcf_headers.append(line)
                continue
            fields = line.split("\t")
            variants[fields[0], int(fields[1])] = {
                "chrom": fields[0], "pos": int(fields[1]),
                "ref": fields[3], "alt": fields[4],
                "genotype": fields[9].split(":")[0], "vcf_record": line,
            }

    manifest = {"schema_version": 1, "datasets": {}}
    for name, filename in SOURCES.items():
        with pysam.AlignmentFile(source_dir / filename) as alignment:
            header = alignment.text
            reads = list(alignment)
        cases = [variants[key] for key in CASES[name]]
        if name == "illumina_star":
            # All observations near the four chosen loci, plus their mates and
            # alternative alignments available in the published chr4 subset.
            names = {r.query_name for r in reads if any(
                r.get_overlap(v["pos"] - 1, v["pos"] + len(v["ref"])) > 0
                for v in cases)}
            # Keep three genuine multimapped templates to pin NH/MAPQ behavior.
            names.update(sorted({r.query_name for r in reads if r.get_tag("NH") > 1})[:3])
            # One directly inspected overlapping pair, outside the benchmark
            # loci, for a fragment-merging regression (not a variant claim).
            names.add("SRR1258218.22918678")
            selected = [r for r in reads if r.query_name in names]
        elif name == "pacbio_minimap2":
            # A small stratified edge-case collection, NOT an abundance sample.
            groups = Counter()
            selected = []
            for read in reads:
                qualities = read.query_qualities
                key = (read.flag, tuple(sorted(set(qualities or []))))
                if groups[key] < 4:
                    selected.append(read)
                    groups[key] += 1
        else:
            # Preserve the entire original region, including low-MAPQ reads,
            # supplementary records and reads spliced across the queried loci.
            selected = reads
        raw = (header + "".join(r.to_string() + "\n" for r in selected)).encode("ascii")
        target = output_dir / (name + ".sam.gz")
        # Explicit mtime and empty filename make gzip output reproducible.
        with target.open("wb") as handle:
            with gzip.GzipFile(filename="", mode="wb", fileobj=handle, mtime=0) as compressed:
                compressed.write(raw)
        manifest["datasets"][name] = {
            "file": target.name,
            "sam_sha256": sha256(raw).hexdigest(),
            "source_region_records": len(reads),
            "fixture_records": len(selected),
            "records": [record_digest(r) for r in selected],
            "mapq_counts": dict(sorted(Counter(r.mapping_quality for r in selected).items())),
            "variants": [dict(v, observations=[obs for r in selected
                if (obs := cigar_observation(r, v)) is not None]) for v in cases],
        }
    (output_dir / "manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
    selected_records = sorted({v["vcf_record"] for key in CASES.values() for v in
                               (variants[k] for k in key)},
                              key=lambda line: (line.split("\t")[0], int(line.split("\t")[1])))
    (output_dir / "HG001.selected.vcf").write_text(
        "\n".join(vcf_headers + selected_records) + "\n")
    selected_bed = []
    with (source_dir / "HG001.bed").open() as bed:
        for line in bed:
            chrom, start, end = line.rstrip().split("\t")[:3]
            if any(chrom == v["chrom"] and int(start) <= v["pos"] - 1
                   and int(end) >= v["pos"] - 1 + len(v["ref"])
                   for group in manifest["datasets"].values() for v in group["variants"]):
                selected_bed.append(line.rstrip())
    (output_dir / "HG001.selected.bed").write_text("\n".join(selected_bed) + "\n")
    print(json.dumps({name: {"reads": data["fixture_records"],
                             "mapq": data["mapq_counts"]}
                      for name, data in manifest["datasets"].items()}, indent=2))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("source_dir", type=Path)
    parser.add_argument("output_dir", type=Path)
    args = parser.parse_args()
    build(args.source_dir, args.output_dir)
