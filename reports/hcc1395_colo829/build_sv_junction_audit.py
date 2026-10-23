"""Replay selected original-read SV candidates without a MAPQ discovery gate.

The 500-bp screen is a triage step. Low MAPQ on a split alignment is placement
uncertainty, not absence of its sequenced junction. This audit separates a
sequence-supported RNA join, a MAPQ sensitivity subset, mate context, and
regional path assembly. It does not infer a coding frame without models.
"""

import argparse
import csv
from pathlib import Path

import pysam

from isovar.fusion import FusionBreakpoint
from isovar.read_collector import ReadCollector
from isovar.sv_rna import reconstruct_sv_rna


HERE = Path(__file__).parent
ROOT = HERE.parents[1]
HCC_BAM = ROOT / "tests/data/hcc1395_sv_junction/SRR23456709.sv_candidates.grch38.bam"
COLO_BAM = ROOT / "tests/data/colo829_whole_run/SRR8615617.all_sites.grch38.bam"
OUTPUT = HERE / "sv_junction_audit.tsv"

# Oriented hypotheses checked against original records. These are a curated
# subset of the screen, not an automatic interpretation of every VCF SV type.
EVENTS = (
    ("HCC1395", "chr13_52688329_chr13_60966095_+-", "DEL",
     ("chr13", 52688329, "+"), ("chr13", 60966095, "+"), 5, 0),
    ("HCC1395", "chr21_42975947_chr21_43035721_-+", "DUP",
     ("chr21", 43035721, "+"), ("chr21", 42975947, "+"), 2, 0),
    ("HCC1395", "chr9_83969153_chr9_84015113_-+", "DUP",
     ("chr9", 83969153, "-"), ("chr9", 84015113, "-"), 7, 0),
    ("COLO829", "COLO829_SOMATIC_SV3_DEL", "DEL",
     ("chr1", 224595102, "+"), ("chr1", 224598323, "+"), 0, 6),
    ("COLO829", "COLO829_SOMATIC_SV42_DEL", "DEL",
     ("chr16", 58590638, "+"), ("chr16", 58629419, "+"), 0, 1),
    ("COLO829", "COLO829_SOMATIC_SV49_DEL", "DEL",
     ("chr22", 33363264, "+"), ("chr22", 33442413, "+"), 0, 1),
)
FIELDS = (
    "sample", "call_id", "class", "rna_bam", "screen_split_candidates",
    "screen_mate_candidates", "junction_fragments_all_mapq",
    "junction_fragments_mapq20_subset", "junction_read_ids_all_mapq",
    "junction_piece_mapqs_all_mapq", "junction_left", "junction_right",
    "junction_unplaced_bases", "assembled_event_paths_retained",
    "max_assembled_context_bases", "assembly_limitations", "frame_status",
)


def _reconstruct(bam, event, min_mapping_quality):
    sample, call_id, _, donor, acceptor, _, _ = event
    with pysam.AlignmentFile(str(bam), "rb") as alignment:
        return reconstruct_sv_rna(
            alignment, event_id=call_id, reference_name="GRCh38",
            sample_id=sample, source=str(bam),
            donor=FusionBreakpoint(*donor), acceptor=FusionBreakpoint(*acceptor),
            regions=[], references=[], event_provenance={"dna_call_id": call_id},
            read_collector=ReadCollector(
                min_mapping_quality=min_mapping_quality, use_soft_clipped_bases=True),
            max_records=5000, max_paths=10, max_queries=100,
            assemble=True, include_regional_candidates=False,
            dense_support=True)


def _breakpoint_paths(result):
    return [(path, junction) for path in result["paths"] for junction in path["junctions"]
            if junction["relation"] == "breakpoint_junction"]


def audit_event(event, bam):
    """Return distinct-fragment junction counts and bounded assembly context."""
    sample, call_id, cls, _, _, split, mates = event
    full = _reconstruct(bam, event, 0)
    subset = _reconstruct(bam, event, 20)
    joins = _breakpoint_paths(full)
    high_mapq_joins = _breakpoint_paths(subset)
    # The same junction can occur on several paths. Never add their support.
    best = max((junction for _, junction in joins),
               key=lambda j: (j["direct_support"]["fragments"], repr(j["left"]), repr(j["right"])),
               default=None)
    observations = full["observations"]
    witnesses = (sorted({observations[key]["identity"][1]
                         for key in best["direct_observations"]}) if best else [])
    record_ids = (set().union(*(set(observations[key]["records"])
                                for key in best["direct_observations"])) if best else set())
    piece_mapqs = sorted({full["record_evidence"][record_id]["mapping_quality_raw"]
                          for record_id in record_ids})
    event_paths = {path["path_id"]: path for path, _ in joins}
    return dict(
        sample=sample, call_id=call_id, **{"class": cls}, rna_bam=Path(bam).name,
        screen_split_candidates=split, screen_mate_candidates=mates,
        junction_fragments_all_mapq=best["direct_support"]["fragments"] if best else 0,
        junction_fragments_mapq20_subset=max(
            (j["direct_support"]["fragments"] for _, j in high_mapq_joins), default=0),
        junction_read_ids_all_mapq=";".join(witnesses),
        junction_piece_mapqs_all_mapq=";".join(map(str, piece_mapqs)),
        junction_left=":".join(map(str, best["left"])) if best else "",
        junction_right=":".join(map(str, best["right"])) if best else "",
        junction_unplaced_bases=best["unplaced_bases"] if best else "",
        assembled_event_paths_retained=len(event_paths),
        max_assembled_context_bases=max((len(path["sequence"]) for path in event_paths.values()), default=0),
        assembly_limitations=";".join(full["limitations"]),
        frame_status="reference_models_unavailable",
    )


def write_audit(output=OUTPUT, hcc_bam=HCC_BAM, colo_bam=COLO_BAM):
    """Write the selected-event table from source-derived indexed BAMs."""
    with Path(output).open("w", newline="") as stream:
        writer = csv.DictWriter(stream, FIELDS, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for event in EVENTS:
            writer.writerow(audit_event(event, hcc_bam if event[0] == "HCC1395" else colo_bam))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--hcc-bam", type=Path, default=HCC_BAM)
    parser.add_argument("--colo-bam", type=Path, default=COLO_BAM)
    parser.add_argument("--output", type=Path, default=OUTPUT)
    args = parser.parse_args()
    write_audit(args.output, args.hcc_bam, args.colo_bam)


if __name__ == "__main__":
    main()
