"""Screen original COLO829 RNA for observed alignments near nominated SV ends.

This is an event-level evidence screen, not structural RNA reconstruction.
Two observed alignments from one read can be a junction candidate, but their
sequence, strand, splice alternative, and ORF still require reconciliation.
"""

import argparse
from collections import defaultdict
import csv
import gzip
from pathlib import Path

import pysam

from .build_atlas_regions import breakpoints


HERE = Path(__file__).parent
FIXTURE = HERE.parents[1] / "tests/data/colo829_whole_run"
SOURCES = ("PaulinPaz_2024_v2", "ValleInclan_2022_lifted")
FIELDS = (
    "sample", "dna_source", "class", "call_id", "chrom1", "pos1", "chrom2",
    "pos2", "vcf_alt", "dna_filter", "dna_support", "rna_run", "window_bp",
    "end1_alignment_records",
    "end2_alignment_records", "end1_fragments", "end2_fragments",
    "observed_split_alignment_candidates", "candidate_read_ids",
    "cross_end_mate_pair_candidates", "pair_candidate_read_ids", "rna_status",
    "path_status", "protein_status",
)


def source_sv_rows(vcfs):
    """Yield original VCF rows without merging the two benchmark definitions."""
    for source_name, path in zip(SOURCES, vcfs):
        with path.open() as source:
            for line in source:
                if line.startswith("#"):
                    continue
                cols = line.rstrip("\n").split("\t")
                coords = breakpoints(cols)
                cls = next((part.split("=", 1)[1] for part in cols[7].split(";")
                            if part.startswith("SVTYPE=")), "SV")
                support = next((part.split("=", 1)[1] for part in cols[7].split(";")
                                if part.startswith("SUPPORT=")), "")
                yield source_name, cls, cols[2], coords, cols[4], cols[6], support


def window_alignments(alignment, chrom, pos, flank):
    """Return high-mapq observed alignments around a one-based breakend."""
    contig = chrom if chrom.startswith("chr") else f"chr{chrom}"
    return [read for read in alignment.fetch(contig, max(0, pos - flank - 1),
                                             pos + flank)
            if read.mapping_quality >= 20 and not read.is_unmapped
            and not read.is_secondary and not read.is_duplicate
            and not read.is_qcfail]


def split_candidates(first, second):
    """Find QNAMEs with two observed segments, including a supplementary one."""
    by_end = []
    for reads in (first, second):
        indexed = defaultdict(list)
        for read in reads:
            indexed[(read.get_tag("RG") if read.has_tag("RG") else "", read.query_name)].append(read)
        by_end.append(indexed)
    candidates = []
    for key in by_end[0].keys() & by_end[1].keys():
        if any((a.is_supplementary or b.is_supplementary) and
               a.is_read1 == b.is_read1 and a.is_read2 == b.is_read2 and
               (a.reference_name, a.reference_start, a.flag) !=
               (b.reference_name, b.reference_start, b.flag)
               for a in by_end[0][key] for b in by_end[1][key]):
            candidates.append(key[1])
    return sorted(set(candidates))


def cross_end_mate_candidates(first, second):
    """Find pair names with opposite primary mates at the two SV windows."""
    first_read1 = {r.query_name for r in first if r.is_read1 and not r.is_supplementary}
    first_read2 = {r.query_name for r in first if r.is_read2 and not r.is_supplementary}
    second_read1 = {r.query_name for r in second if r.is_read1 and not r.is_supplementary}
    second_read2 = {r.query_name for r in second if r.is_read2 and not r.is_supplementary}
    return sorted((first_read1 & second_read2) | (first_read2 & second_read1))


def screen_rows(vcfs, alignment, flank=500):
    """Yield one explicit RNA status per source VCF record."""
    for source, cls, call_id, coords, vcf_alt, dna_filter, dna_support in source_sv_rows(vcfs):
        chrom1, pos1 = coords[0]
        chrom2, pos2 = coords[1] if len(coords) > 1 else ("", 0)
        first = window_alignments(alignment, chrom1, pos1, flank)
        second = window_alignments(alignment, chrom2, pos2, flank) if chrom2 else []
        witnesses = split_candidates(first, second) if chrom2 else []
        distinct_windows = chrom1 != chrom2 or abs(pos1 - pos2) > 2 * flank
        pair_witnesses = (cross_end_mate_candidates(first, second)
                          if chrom2 and distinct_windows else [])
        if not chrom2:
            status = "unresolved_partner_coordinate"
        elif witnesses:
            status = "split_alignment_candidate_unresolved"
        elif pair_witnesses:
            status = "cross_end_pair_context_unresolved"
        elif first and second:
            status = "both_ends_covered_no_split_candidate"
        elif first or second:
            status = "one_end_covered_no_split_candidate"
        else:
            status = "no_breakend_window_alignment"
        yield {
            "sample": "COLO829", "dna_source": source, "class": cls,
            "call_id": call_id, "chrom1": chrom1, "pos1": pos1,
            "chrom2": chrom2, "pos2": pos2 or "NA",
            "vcf_alt": vcf_alt,
            "dna_filter": dna_filter, "dna_support": dna_support,
            "rna_run": "SRR8615617",
            "window_bp": flank, "end1_alignment_records": len(first),
            "end2_alignment_records": len(second),
            "end1_fragments": len({r.query_name for r in first}),
            "end2_fragments": len({r.query_name for r in second}),
            "observed_split_alignment_candidates": len(witnesses),
            "candidate_read_ids": ";".join(witnesses), "rna_status": status,
            "cross_end_mate_pair_candidates": len(pair_witnesses),
            "pair_candidate_read_ids": ";".join(pair_witnesses),
            "path_status": "unresolved", "protein_status": "not_called_unresolved_path",
        }


def write_screen(vcfs, bam, output, inventory=HERE / "somatic_call_inventory.tsv.gz",
                 flank=500):
    """Write and count the source-separated SV screen rows."""
    with gzip.open(inventory, "rt", newline="") as source:
        inventory_keys = {(row["source"], row["class"], row["call_id"])
                          for row in csv.DictReader(source, delimiter="\t")
                          if row["sample"] == "COLO829" and row["source"] in SOURCES}
    seen = set()
    with pysam.AlignmentFile(str(bam), "rb") as alignment, output.open("w", newline="") as target:
        writer = csv.DictWriter(target, FIELDS, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        count = 0
        for row in screen_rows(vcfs, alignment, flank):
            key = (row["dna_source"], row["class"], row["call_id"])
            if key not in inventory_keys or key in seen:
                raise ValueError(f"SV source record does not uniquely match inventory: {key}")
            seen.add(key)
            writer.writerow(row)
            count += 1
    if seen != inventory_keys:
        raise ValueError(f"SV inventory has {len(inventory_keys - seen)} unreported records")
    return count


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--sv-vcf", type=Path, nargs=2,
                        default=(FIXTURE / "colo829-sv-grch38.vcf",
                                 FIXTURE / "colo829-sv-2022-hg38lifted.vcf"))
    parser.add_argument("--inventory", type=Path, default=HERE / "somatic_call_inventory.tsv.gz")
    parser.add_argument("--bam", type=Path,
                        default=FIXTURE / "SRR8615617.all_sites.grch38.bam")
    parser.add_argument("--output", type=Path, default=HERE / "colo829_sv_rna_screen.tsv")
    args = parser.parse_args()
    print(f"Screened {write_screen(args.sv_vcf, args.bam, args.output, args.inventory)} SV records")


if __name__ == "__main__":
    main()
