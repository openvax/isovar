"""Annotate nominated HCC1395/COLO829 SV footprints against Ensembl 87 exons.

Intrachromosomal DEL/DUP/INV records query the full nominated interval. An
insertion, translocation or breakend queries its nominated breakpoint positions.
This is genomic involvement, not proof of an expressed exon junction.
"""

import csv
import gzip
import hashlib
import io
import json
from collections import Counter
from pathlib import Path

from pyensembl import cached_release

from isovar.evidence_report import read_tsv


HERE = Path(__file__).parent
OUTPUT = HERE / "sv_exon_annotations.tsv.gz"
METADATA = HERE / "sv_exon_annotations.tsv.metadata.json"
FIELDS = ("sample", "dna_source", "class", "call_id", "chrom1", "pos1",
          "chrom2", "pos2", "exon_overlap", "exon_rule")
INPUTS = (HERE / "hcc1395_direct_rna_sv_screen.tsv",
          HERE / "colo829_sv_rna_screen.tsv")


def sv_footprints(row):
    """Return 1-based inclusive intervals and their declared footprint rule."""
    first = row["chrom1"].removeprefix("chr"), int(row["pos1"])
    second = ((row["chrom2"].removeprefix("chr"), int(row["pos2"]))
              if row["pos2"].isdigit() and row["chrom2"] not in {"", "NA"} else None)
    if row["class"] == "INS":
        positions = {first, (first[0], first[1] + 1)}
        if second is not None:
            positions.add(second)
        return [(chrom, pos, pos) for chrom, pos in sorted(positions)], "insertion_flanks"
    if second is None:
        return [(first[0], first[1], first[1])], "one_nominated_breakpoint_only"
    if row["class"] in {"DEL", "DUP", "INV"} and first[0] == second[0]:
        return [(first[0], min(first[1], second[1]),
                 max(first[1], second[1]))], "intrachromosomal_span"
    positions = sorted({first, second})
    return [(chrom, pos, pos) for chrom, pos in positions], "nominated_breakpoints"


def exon_hit(genome, intervals):
    """Return whether any Ensembl exon intersects one of the SV footprints."""
    return any(bool(genome.db.run_sql_query(
        'SELECT 1 FROM exon WHERE seqname = ? AND start <= ? AND end >= ? LIMIT 1',
        query_params=[chrom, end, start])) for chrom, start, end in intervals)


def build(output=OUTPUT, genome=None):
    genome = genome or cached_release(87)
    rows = []
    for source in INPUTS:
        for row in read_tsv(source):
            intervals, rule = sv_footprints(row)
            hit = exon_hit(genome, intervals)
            overlap = "True" if hit else "Unknown" if rule == "one_nominated_breakpoint_only" else "False"
            rows.append({field: row[field] for field in FIELDS[:8]} |
                        dict(exon_overlap=overlap, exon_rule=rule))
    rows.sort(key=lambda row: (row["sample"], row["dna_source"], row["call_id"]))
    output.parent.mkdir(parents=True, exist_ok=True)
    with output.open("wb") as raw, gzip.GzipFile(filename="", fileobj=raw,
                                                  mode="wb", mtime=0) as zipped:
        with io.TextIOWrapper(zipped, encoding="utf-8", newline="") as target:
            writer = csv.DictWriter(target, FIELDS, delimiter="\t", lineterminator="\n")
            writer.writeheader()
            writer.writerows(rows)
    if output == OUTPUT:
        metadata = {
            "schema": "isovar.sv_exon_annotations.v1",
            "genome": "human Ensembl 87 GRCh38",
            "coordinate_system": "1-based inclusive",
            "source_sha256": {path.name: hashlib.sha256(path.read_bytes()).hexdigest()
                              for path in INPUTS},
            "records_by_sample_and_overlap": {
                sample: dict(sorted(Counter(row["exon_overlap"] for row in rows
                                            if row["sample"] == sample).items()))
                for sample in ("HCC1395", "COLO829")},
        }
        METADATA.write_text(json.dumps(metadata, indent=2, sort_keys=True) + "\n")
    return len(rows)


if __name__ == "__main__":
    print("Annotated SV source records:", build())
