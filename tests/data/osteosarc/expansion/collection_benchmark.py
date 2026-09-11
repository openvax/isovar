"""Full-original-region #229 collection benchmark and allocation attribution.

Run in a fresh process. Allocation-instrumented runs are NOT timing/RSS
benchmarks: tracer and snapshot overhead are reported explicitly. End-to-end
RNA/protein benchmarks use the separate existing ``benchmark`` command.
"""

import argparse
from hashlib import sha256
import json
import logging
import math
from pathlib import Path
import platform
import sys
import time
import tracemalloc

import pysam
from varcode import Variant

import isovar
from isovar.allele_read_helpers import allele_reads_from_locus_reads
from isovar.read_collector import ReadCollector
from isovar.read_evidence import ReadEvidence
from isovar.variant_helpers import base0_interval_for_variant_fields, trim_variant
from tests.data.osteosarc.expansion.benchmark import peak_rss_bytes
from tests.data.osteosarc.expansion.inventory import digest, write_json
from tests.data.osteosarc.expansion.runner import counts_from_evidence, primary_alignment, time_limit

# Patch this module-level name in tests: rebinding tracemalloc.take_snapshot on
# the stdlib module would replace it for every other test in the process.
take_snapshot = tracemalloc.take_snapshot


def read_fingerprint(reads):
    """Hash every field in every original object, preserving order/multiplicity.

    Sequence-valued fields also record their container type: JSON writes a list,
    a tuple and an array identically, which would hide a representation change.
    """
    checksum = sha256()
    for read in reads:
        values = {}
        for name in read._fields:
            value = getattr(read, name)
            if not isinstance(value, (str, int, float, type(None))):
                values[name + "_container"] = type(value).__name__
                value = list(value)
            values[name] = value
        checksum.update(json.dumps(values, sort_keys=True, separators=(",", ":")).encode() + b"\n")
    return checksum.hexdigest()


def run(args):
    if not math.isfinite(args.timeout) or not 0 < args.timeout <= 1800:
        raise ValueError("Use a finite timeout greater than zero and at most 1800 seconds")
    if tracemalloc.is_tracing():
        raise ValueError("Run this benchmark in a fresh process without an existing allocation tracer")
    started = time.perf_counter()
    output = Path(args.output)
    output.mkdir(parents=True, exist_ok=False)
    destination = Path(args.destination)
    source = destination / "alignments" / args.source_id
    bam = source / "regions-GRCh38.bam"
    receipt = source / "regions-GRCh38.json"
    original_sha = digest(bam)
    if original_sha != json.loads(receipt.read_text())["bam_sha256"]:
        raise ValueError("Regional BAM checksum mismatch")
    input_sha = original_sha  # Hash each multi-GB regional BAM only once.
    if args.mode == "primary_only":
        bam = primary_alignment(bam, original_sha)
        input_sha = digest(bam)
    inventory = destination / "inventory-GRCh38-validated.json"
    record = next(r for r in json.loads(inventory.read_text())["variants"] if r["variant_id"] == args.variant_id)
    variant = Variant(record["chrom"], record["pos"], record["ref"], record["alt"])
    position, ref, alt = trim_variant(variant)
    base0_start, base0_end = base0_interval_for_variant_fields(position, ref, alt)
    root = Path(isovar.__file__).resolve().parent.parent
    identity = dict(
        source_id=args.source_id, variant_id=args.variant_id, mode=args.mode,
        source_bam_sha256=original_sha, input_bam_sha256=input_sha,
        input_index_sha256=digest(str(bam) + ".bai"), receipt_sha256=digest(receipt),
        inventory_sha256=digest(inventory), benchmark_sha256=digest(__file__),
        source_files={str(p.relative_to(root)): digest(p) for p in sorted((root / "isovar").glob("*.py"))},
        isovar=isovar.__version__, python=platform.python_version(), platform=platform.platform(),
        allocation_instrumented=args.allocations, merge_overlapping_fragments=True,
        timeout_seconds=args.timeout, setup_seconds=time.perf_counter() - started,
    )
    write_json(output / "identity.json", identity)

    class ProfiledCollector(ReadCollector):
        @classmethod
        def _merge_overlapping_locus_reads(cls, reads):
            if args.allocations:
                current, peak = tracemalloc.get_traced_memory()
                overhead = tracemalloc.get_tracemalloc_memory()
                snapshot = take_snapshot()
                # Stop tracking before aggregating millions of trace entries:
                # the analysis itself is not a collection allocation, and
                # keeping tracer bookkeeping alive can exhaust host memory.
                tracemalloc.stop()
                stats = snapshot.statistics("lineno")
                allocations = dict(
                    stage="before_mate_merging", locus_read_count=len(reads),
                    query_bases=sum(len(r.sequence) for r in reads),
                    traced_current_bytes=current, traced_peak_bytes=peak,
                    tracer_overhead_bytes=overhead, process_peak_rss_bytes=peak_rss_bytes(),
                    locations=[dict(location=str(s.traceback), bytes=s.size, count=s.count) for s in stats[:25]])
                write_json(output / "allocations.json", allocations)
                print(json.dumps(allocations, sort_keys=True), flush=True)
            return super()._merge_overlapping_locus_reads(reads)

    disabled = logging.root.manager.disable
    logging.disable(logging.CRITICAL)
    result = dict(status="running")
    collection_started, cpu_started = time.perf_counter(), time.process_time()
    try:
        if args.allocations:
            tracemalloc.start(1)
        # One budget bounds collection and the post-collection work after it;
        # fingerprinting and conversion are themselves unbounded on deep loci.
        with time_limit(args.timeout):
            with pysam.AlignmentFile(bam) as handle:
                chromosome = ReadCollector._infer_chromosome_name(record["chrom"], handle.references)
                if chromosome is None:
                    raise ValueError("Variant contig not present in source BAM")
                reads = ProfiledCollector(merge_overlapping_fragments=True).get_locus_reads(
                    handle, chromosome, base0_start, base0_end, position, ref, alt)
            completed = dict(
                locus_collection_wall_seconds=time.perf_counter() - collection_started,
                locus_collection_cpu_seconds=time.process_time() - cpu_started,
                locus_collection_peak_rss_bytes=peak_rss_bytes(), locus_read_count=len(reads))
            # Fingerprinting and conversion are explicit post-collection work,
            # not part of the collection measurement. The whole-process peak
            # follows.
            post_started = time.perf_counter()
            completed["locus_reads_sha256"] = read_fingerprint(reads)
            allele_reads = allele_reads_from_locus_reads(reads)
            completed["allele_reads_sha256"] = read_fingerprint(allele_reads)
            completed["allele_read_count"] = len(allele_reads)
            completed["counts"] = counts_from_evidence(
                ReadEvidence.from_variant_and_allele_reads(variant, allele_reads))
            completed.update(post_collection_wall_seconds=time.perf_counter() - post_started,
                             whole_process_peak_rss_bytes=peak_rss_bytes())
        # Record completed evidence in one step: an interruption at any earlier
        # point leaves an error record with no counts, never partial evidence.
        result.update(completed, status="ok")
    except BaseException as error:
        result.update(status="error", error=f"{type(error).__name__}: {error}")
        raise
    finally:
        if args.allocations and tracemalloc.is_tracing():
            tracemalloc.stop()
        logging.disable(disabled)
        try:
            write_json(output / "result.json", result)
        except OSError:
            # A full disk must not hide the run's own error or lose its record:
            # keep the record in the log, and re-raise only if nothing else is.
            print(json.dumps(result, sort_keys=True), file=sys.stderr, flush=True)
            if result["status"] != "error":
                raise
    print(json.dumps(result, sort_keys=True), flush=True)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--destination", required=True)
    parser.add_argument("--source-id", required=True)
    parser.add_argument("--variant-id", default="MT_ND5-chrM-12994")
    parser.add_argument("--mode", choices=["defaults", "primary_only"], default="defaults")
    parser.add_argument("--output", required=True)
    parser.add_argument("--timeout", type=float, default=600)
    parser.add_argument("--allocations", action="store_true")
    run(parser.parse_args())
