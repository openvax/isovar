"""Time exact dense recount on a COLO829 SV3 candidate BAM."""

import argparse
from contextlib import closing
import gzip
import json
from pathlib import Path
import sqlite3
from tempfile import TemporaryDirectory
from time import monotonic

import pysam

from isovar.fusion import FusionBreakpoint, FusionReference
from isovar.read_collector import ReadCollector
from isovar.sv_rna import reconstruct_sv_rna
from isovar.sv_rna_dense import _DiskGroups, _spool
from reports.hcc1395_colo829.build_sv_junction_audit import EVENTS


ROOT = Path(__file__).resolve().parents[2]
DATA = ROOT / "tests/data/colo829_whole_run"


def benchmark_templates(bam_path, scratch_dir=None):
    """Time source indexing and one complete paired-template traversal."""
    with TemporaryDirectory(prefix="isovar-sv-benchmark-", dir=scratch_dir) as directory:
        with pysam.AlignmentFile(str(bam_path), "rb") as alignment, closing(
                sqlite3.connect(str(Path(directory) / "records.sqlite"))) as db:
            indexed = monotonic()
            scanned, excluded, unique, segments, source = _spool(
                alignment, db, ReadCollector(min_mapping_quality=0, use_soft_clipped_bases=True),
                directory)
            spool_seconds = round(monotonic() - indexed, 2)
            traversed = monotonic()
            with pysam.AlignmentFile(source, "rb") as reader:
                groups = _DiskGroups(db, reader, alignment.header, segments)
                templates = records = 0
                for template in groups.templates():
                    templates += 1
                    records += len(template)
            return dict(
                spool_seconds=spool_seconds,
                traversal_seconds=round(monotonic() - traversed, 2),
                records_scanned=scanned, unique_records=unique,
                excluded_records=excluded, eligible_segments=segments,
                templates=templates, eligible_records=records)


def benchmark(bam_path, scratch_dir=None, max_records=20000):
    """Return wall time and exact join support for one original-read replay."""
    models = json.loads(gzip.open(DATA / "ensembl87_cnih3_models.json.gz", "rt").read())
    references = [FusionReference(**row) for row in models["models"]]
    sample, call_id, _, donor, acceptor, *_ = EVENTS[3]
    started = monotonic()
    with pysam.AlignmentFile(str(bam_path), "rb") as alignment:
        result = reconstruct_sv_rna(
            alignment, event_id=call_id, reference_name="GRCh38", sample_id=sample,
            source=str(bam_path), donor=FusionBreakpoint(*donor),
            acceptor=FusionBreakpoint(*acceptor), regions=[], references=references,
            event_provenance={"dna_call_id": call_id},
            read_collector=ReadCollector(min_mapping_quality=0, use_soft_clipped_bases=True),
            max_records=max_records, max_queries=200, max_paths=20,
            max_extension_segments=1200, assemble=True,
            include_regional_candidates=False, dense_support=True,
            exon_anchored_clips=True, scratch_dir=scratch_dir)
    joined = [(path, junction) for path in result["paths"]
              for junction in path["junctions"]
              if junction["relation"] == "breakpoint_junction"]
    return dict(
        wall_seconds=round(monotonic() - started, 2),
        records_scanned=result["support_acquisition"]["records_scanned"],
        eligible_segments=result["support_acquisition"]["eligible_segments"],
        selected_records=result["discovery"]["selected_records"],
        path_count=len(result["paths"]),
        joined_paths=[dict(path_id=path["path_id"],
                           direct_fragments=junction["direct_support"]["fragments"],
                           local_exact_fragments=path.get("event_local_support", {}).get("fragments"))
                      for path, junction in joined],
        limitations=result["limitations"])


def main():
    """Run the pinned replay, or benchmark a larger BAM with the same call."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--bam", type=Path,
                        default=DATA / "SRR8615617.cnih3_source_candidates.grch38.bam")
    parser.add_argument("--scratch-dir", type=Path)
    parser.add_argument("--max-records", type=int, default=20000)
    parser.add_argument("--templates-only", action="store_true")
    args = parser.parse_args()
    result = (benchmark_templates(args.bam, args.scratch_dir) if args.templates_only
              else benchmark(args.bam, args.scratch_dir, args.max_records))
    print(json.dumps(result, sort_keys=True))


if __name__ == "__main__":
    main()
