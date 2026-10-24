"""Benchmark bounded SV discovery against pinned original RNA alignments.

Run one case per process so peak RSS is attributable to that case. The input
may be a larger local original-read subset; no reads are synthesized here.
"""

import argparse
import gzip
from hashlib import sha256
import json
from pathlib import Path
import resource
import sys
from tempfile import TemporaryDirectory
from threading import Event, Thread
from time import monotonic

import pysam

from isovar.fusion import FusionBreakpoint, FusionReference
from isovar.read_collector import ReadCollector
from isovar.sv_rna import reconstruct_sv_rna
from reports.hcc1395_colo829.build_sv_junction_audit import EVENTS


ROOT = Path(__file__).resolve().parents[2]
SOURCES = {
    "colo829": (ROOT / "tests/data/colo829_whole_run/SRR8615617.cnih3_source_candidates.grch38.bam",
                ROOT / "tests/data/colo829_whole_run/ensembl87_cnih3_models.json.gz", 3, "models"),
    "hcc1395": (ROOT / "tests/data/hcc1395_sv_junction/SRR23456709.sv_candidates.grch38.bam",
                ROOT / "tests/data/hcc1395_sv_junction/ensembl87_models.json.gz", 0, "SUGT1"),
}


def _scratch_bytes(root):
    size = 0
    for path in root.rglob("*"):
        try:
            size += path.stat().st_size if path.is_file() else 0
        except FileNotFoundError:
            pass  # Reconstruction may remove a temporary file during sampling.
    return size


def benchmark(sample, bam_path, max_records, max_evidence_groups, all_aligned_regions=False):
    """Measure one complete discovery and original-input support recount."""
    _, model_path, event_index, model_key = SOURCES[sample]
    models = json.loads(gzip.open(model_path, "rt").read())["models"]
    if model_key != "models":
        models = models[model_key]
    references = [FusionReference(**row) for row in models]
    sample_name, call_id, _, donor, acceptor, *_ = EVENTS[event_index]
    stop = Event()
    peak_scratch = [0]
    with TemporaryDirectory(prefix="isovar-group-benchmark-") as directory:
        scratch = Path(directory)

        def watch_scratch():
            while not stop.wait(0.1):
                peak_scratch[0] = max(peak_scratch[0], _scratch_bytes(scratch))

        watcher = Thread(target=watch_scratch, daemon=True)
        watcher.start()
        started = monotonic()
        try:
            with pysam.AlignmentFile(str(bam_path), "rb") as alignment:
                regions = ([(name, 0, length) for name, length in zip(
                    alignment.references, alignment.lengths) if length]
                           if all_aligned_regions else [])
                result = reconstruct_sv_rna(
                    alignment, event_id=call_id, reference_name="GRCh38", sample_id=sample_name,
                    source=str(bam_path), donor=FusionBreakpoint(*donor),
                    acceptor=FusionBreakpoint(*acceptor), regions=regions, references=references,
                    event_provenance={"dna_call_id": call_id},
                    read_collector=ReadCollector(min_mapping_quality=0, use_soft_clipped_bases=True),
                    max_records=max_records, max_evidence_groups=max_evidence_groups,
                    max_queries=200, max_paths=20, max_extension_segments=1200,
                    assemble=True, include_regional_candidates=False, dense_support=True,
                    exon_anchored_clips=True, scratch_dir=directory)
        finally:
            stop.set()
            watcher.join()
        elapsed = monotonic() - started
    paths = [dict(path_id=path["path_id"], sequence_sha256=sha256(path["sequence"].encode()).hexdigest(),
                  length=len(path["sequence"]), frame_status=path["frame_status"],
                  full_interval_fragments=path["full_interval_support"]["fragments"],
                  direct_junction_fragments=max((j["direct_support"]["fragments"] for j in path["junctions"]),
                                                default=0))
             for path in result["paths"]]
    peak_rss = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    return dict(sample=sample, source=str(bam_path), max_records=max_records,
                all_aligned_regions=all_aligned_regions,
                max_evidence_groups=max_evidence_groups,
                wall_seconds=round(elapsed, 2),
                peak_rss_bytes=peak_rss if sys.platform == "darwin" else peak_rss * 1024,
                peak_scratch_bytes=peak_scratch[0], discovery=result["discovery"],
                support_acquisition=result["support_acquisition"],
                limitations=result["limitations"], paths=paths)


def main():
    """Print a JSON measurement for one source and one pair of limits."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("sample", choices=sorted(SOURCES))
    parser.add_argument("--bam", type=Path)
    parser.add_argument("--max-records", type=int, required=True)
    parser.add_argument("--max-evidence-groups", type=int)
    parser.add_argument("--all-aligned-regions", action="store_true",
                        help="Stress discovery with all source alignments eligible as context")
    args = parser.parse_args()
    bam_path = args.bam or SOURCES[args.sample][0]
    if not bam_path.is_file():
        parser.error("RNA BAM not found: %s" % bam_path)
    print(json.dumps(benchmark(args.sample, bam_path, args.max_records,
                               args.max_evidence_groups, args.all_aligned_regions),
                     sort_keys=True))


if __name__ == "__main__":
    main()
