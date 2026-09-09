"""Limited original-read MT origin diagnostics, never NUMT exclusion or CCF.

The published polymorphic numtA interval overlaps this ND5 locus. A splice
skip must not be mistaken for a continuous read extending beyond that NUMT.
Competitive nuclear/mitochondrial alignment and matched DNA are not performed.
"""

import argparse
from collections import Counter
import json
from pathlib import Path
import sys

import pysam

ROOT = Path(__file__).resolve().parents[4]
sys.path.insert(0, str(ROOT))
from tests.data.osteosarc.expansion.inventory import digest, write_json  # noqa: E402
from tests.data.osteosarc.expansion.runner import PRIMARY_EXCLUDE_FLAGS  # noqa: E402
from tests.real_rna_helpers import cigar_observation  # noqa: E402


NUMTA = (12360, 13227)  # Published MT 12361..13227, converted to half-open.
SOURCES = ("f30f618fb76a0e49", "7010dd389c50467d", "0066232879babe83")


def variant_block(read, position, max_deletion=50):
    """Return the aligned block containing a 0-based variant; N always breaks."""
    blocks, start, end = [], None, None
    pos = read.reference_start
    for operation, length in read.cigartuples or []:
        if operation in (0, 7, 8):
            start = pos if start is None else start
            end = pos + length
        elif operation in (3, 4, 5) or (operation == 2 and length > max_deletion):
            if start is not None:
                blocks.append((start, end))
            start, end = None, None
        if operation in (0, 2, 3, 7, 8):
            pos += length
    if start is not None:
        blocks.append((start, end))
    return next((b for b in blocks if b[0] <= position < b[1]), None)


def source_diagnostic(destination, sid, record):
    directory = Path(destination) / "alignments" / sid
    path = directory / "regions-GRCh38.bam"
    receipt = json.loads(path.with_suffix(".json").read_text())
    if receipt["status"] != "ok" or digest(path) != receipt["bam_sha256"]:
        raise ValueError("MT diagnostic source integrity failure")
    groups, mapqs, lengths, nh = Counter(), Counter(), Counter(), Counter()
    flags, quality = Counter(), Counter()
    blocks = []
    with pysam.AlignmentFile(path) as bam:
        for read in bam.fetch(record["chrom"], record["pos"] - 2, record["pos"] + 1):
            if read.flag & PRIMARY_EXCLUDE_FLAGS or read.mapping_quality < 1:
                continue
            observed = cigar_observation(read, record)
            if observed is None:
                continue
            group = "ref" if observed["allele"] == record["ref"] else "alt" if observed["allele"] == record["alt"] else "other"
            groups[group] += 1
            if group != "alt":
                continue
            mapqs[str(read.mapping_quality)] += 1
            nh[str(read.get_tag("NH")) if read.has_tag("NH") else "absent"] += 1
            lengths[read.query_alignment_length] += 1
            flags["with_SA"] += read.has_tag("SA")
            flags["with_CIGAR_N"] += any(op == 3 for op, _ in read.cigartuples)
            if observed["qualities"] is None:
                quality["missing"] += 1
            else:
                for cutoff in (0, 10, 20, 30):
                    quality[str(cutoff)] += min(observed["qualities"]) >= cutoff
            block = variant_block(read, record["pos"] - 1)
            if block:
                blocks.append(block)
                flags["variant_block_extends_outside_published_numtA"] += block[0] < NUMTA[0] or block[1] > NUMTA[1]
    return dict(source_id=sid, source_bam_sha256=receipt["bam_sha256"],
                independent_exact_primary_alignment_counts={g: groups[g] for g in ("ref", "alt", "other")},
                alt_mapq_histogram=dict(mapqs), alt_NH_histogram=dict(nh), alt_quality_retention=dict(quality),
                alt_aligned_query_length_range=[min(lengths), max(lengths)] if lengths else None,
                alt_variant_block_length_range=[min(b - a for a, b in blocks), max(b - a for a, b in blocks)] if blocks else None,
                alt_diagnostics=dict(flags), numt_origin_status="not_ruled_out", dna_heteroplasmy=None, cancer_cell_fraction=None)


def build(destination, output):
    destination = Path(destination)
    inventory = json.loads((destination / "inventory-GRCh38-validated.json").read_text())
    record = next(v for v in inventory["variants"] if v["chrom"] == "chrM")
    write_json(output, dict(
        variant=record, sources=[source_diagnostic(destination, sid, record) for sid in SOURCES],
        scope="Three diagnostic RNA products, not a NUMT-screened or patient-clonality truth set",
        published_numtA=dict(mt_interval_1based=[12361, 13227], nuclear_insertion="chr21:9676568",
                             source="https://pmc.ncbi.nlm.nih.gov/articles/PMC8896463/"),
        block_definition="N, clipping and deletions over 50 bp break blocks; only the block containing the variant is compared",
        limitations=["Extending outside one published NUMT interval does not exclude other/patient-specific NUMTs or alignment error",
                     "MAPQ 255 is unavailable, not an ordinary Phred probability; NH=1 is not origin proof",
                     "SA presence is only a diagnostic, not a chimeric-origin verdict",
                     "No competitive whole-genome/mitochondrial/NUMT-decoy realignment was performed",
                     "RNA allele fraction is not DNA heteroplasmy or malignant-cell fraction"],
    ))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("destination", type=Path)
    parser.add_argument("output", type=Path)
    args = parser.parse_args()
    build(args.destination, args.output)
