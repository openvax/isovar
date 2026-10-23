"""Replay selected original-read SV candidates with Isovar junction validation.

The 500-bp screen is a triage step. This audit separates an observed,
orientation-compatible RNA join from mate context and from an entire RNA path
seen in one original read. It does not infer a coding frame without models.
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
    "screen_mate_candidates", "mapq1_junction_fragments",
    "mapq20_junction_fragments", "mapq20_junction_read_ids",
    "mapq20_left", "mapq20_right", "mapq20_unplaced_bases",
    "mapq20_event_paths_retained", "mapq20_max_full_path_fragments",
    "mapq20_limitations", "frame_status",
)


def _reconstruct(bam, event, min_mapping_quality, dense_support):
    sample, call_id, _, donor, acceptor, _, _ = event
    with pysam.AlignmentFile(str(bam), "rb") as alignment:
        return reconstruct_sv_rna(
            alignment, event_id=call_id, reference_name="GRCh38",
            sample_id=sample, source=str(bam),
            donor=FusionBreakpoint(*donor), acceptor=FusionBreakpoint(*acceptor),
            regions=[], references=[], event_provenance={"dna_call_id": call_id},
            read_collector=ReadCollector(min_mapping_quality=min_mapping_quality),
            max_records=5000, max_paths=10, max_queries=100,
            assemble=False, include_regional_candidates=False,
            dense_support=dense_support)


def _breakpoint_paths(result):
    return [(path, junction) for path in result["paths"] for junction in path["junctions"]
            if junction["relation"] == "breakpoint_junction"]


def audit_event(event, bam):
    """Return one record; strict counts use all eligible original BAM records."""
    sample, call_id, cls, _, _, split, mates = event
    permissive = _reconstruct(bam, event, 1, False)
    strict = _reconstruct(bam, event, 20, True)
    loose_joins, joins = _breakpoint_paths(permissive), _breakpoint_paths(strict)
    # The same junction can occur on several paths. Never add their support.
    best = max((junction for _, junction in joins),
               key=lambda j: (j["direct_support"]["fragments"], repr(j["left"]), repr(j["right"])),
               default=None)
    observations = strict["observations"]
    witnesses = (sorted({observations[key]["identity"][1]
                         for key in best["direct_observations"]}) if best else [])
    event_paths = {path["path_id"]: path for path, _ in joins}
    return dict(
        sample=sample, call_id=call_id, **{"class": cls}, rna_bam=Path(bam).name,
        screen_split_candidates=split, screen_mate_candidates=mates,
        mapq1_junction_fragments=max((j["direct_support"]["fragments"] for _, j in loose_joins), default=0),
        mapq20_junction_fragments=best["direct_support"]["fragments"] if best else 0,
        mapq20_junction_read_ids=";".join(witnesses),
        mapq20_left=":".join(map(str, best["left"])) if best else "",
        mapq20_right=":".join(map(str, best["right"])) if best else "",
        mapq20_unplaced_bases=best["unplaced_bases"] if best else "",
        mapq20_event_paths_retained=len(event_paths),
        mapq20_max_full_path_fragments=max(
            (path["full_interval_support"]["fragments"] for path in event_paths.values()), default=0),
        mapq20_limitations=";".join(strict["limitations"]),
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
