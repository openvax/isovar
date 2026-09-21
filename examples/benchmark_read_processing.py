"""Offline collection benchmark using only the packaged osteosarc read subset.

Run: python -m examples.benchmark_read_processing output.json
Times exclude acquisition, SAM decoding, annotation and protein reconstruction.
The identical script can run against the parent checkout for comparison.
"""

import argparse
from collections import defaultdict
from hashlib import sha256
import importlib.util
import json
import logging
from pathlib import Path
import platform
import statistics
import time
import tracemalloc
from types import SimpleNamespace

import pysam

from isovar import __version__
from isovar.read_collector import ReadCollector
from isovar.sid_data import BUNDLE, RECIPE, read_json, verify
from isovar.sv_rna_orfs import record_evidence
from isovar.variant_helpers import trim_variant


class Records:
    def __init__(self, reads):
        self.reads = reads

    def fetch(self, chromosome, start, end):
        return (r for r in self.reads if r.reference_name == chromosome
                and r.reference_start < end and r.reference_end > start)


def measure(functions, repeats):
    results = {name: function() for name, function in functions.items()}
    times = {name: dict(seconds=[], cpu_seconds=[]) for name in functions}
    for repeat in range(repeats):
        # Alternate order to reduce drift from load or thermal changes.
        for name in list(functions)[::1 if repeat % 2 else -1]:
            wall, cpu = time.perf_counter(), time.process_time()
            functions[name]()
            times[name]["cpu_seconds"].append(time.process_time() - cpu)
            times[name]["seconds"].append(time.perf_counter() - wall)
    for name, function in functions.items():
        tracemalloc.start()
        function()
        _, peak = tracemalloc.get_traced_memory()
        tracemalloc.stop()
        times[name].update(median_seconds=statistics.median(times[name]["seconds"]),
                           median_cpu_seconds=statistics.median(times[name]["cpu_seconds"]),
                           traced_peak_bytes=peak)
    return results, times


def baseline_module(checkout, module):
    name = "isovar._benchmark_baseline_" + module
    spec = importlib.util.spec_from_file_location(name, checkout / "isovar" / (module + ".py"))
    result = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(result)
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("output", type=Path)
    parser.add_argument("--repeats", type=int, default=7)
    parser.add_argument("--baseline", type=Path, help="Parent checkout for alternating comparisons")
    args = parser.parse_args()
    if args.repeats < 1:
        parser.error("--repeats must be positive")
    logging.disable(logging.CRITICAL)
    root = Path(__file__).resolve().parents[1]
    recipe, bundle = read_json(RECIPE), read_json(BUNDLE)
    verify(recipe, bundle)
    cases = json.loads((root / "tests/data/osteosarc/expansion/corpus/manifest.json").read_text())["cases"]
    inputs = [("osteosarc/expansion/corpus/" + c["bam"], c["variant"]) for c in cases]
    inputs += [("osteosarc/figure_comparisons/corpus/pacbio-indels.bam", c["variant"])
               for c in cases if c["variant"]["pos"] in (3856149, 55627965, 80327830)]
    groups = defaultdict(list)
    for fixture, variant in inputs:
        selection = recipe["fixtures"][fixture]
        source = bundle["sources"][selection["source"]]
        header = pysam.AlignmentHeader.from_dict(source["header"])
        reads = [pysam.AlignedSegment.fromstring(source["records"][key], header)
                 for key in selection["records"]]
        # Platform is a benchmark grouping only, never a collection policy.
        url = source["asset"]["url"].lower()
        technology = "PacBio" if "pacbio" in url else "ONT" if "ont" in url else "Illumina"
        v = SimpleNamespace(contig=variant["chrom"], start=variant["pos"],
                            ref=variant["ref"], alt=variant["alt"])
        pos, ref, alt = trim_variant(v)
        start = pos - 1 if ref else pos
        groups[technology].append((Records(reads), dict(chromosome=variant["chrom"],
            base0_start_inclusive=start, base0_end_exclusive=start + len(ref),
            trimmed_base1_start=pos, trimmed_ref=ref, trimmed_alt=alt)))
    collectors, extractors = dict(current=ReadCollector()), dict(current=record_evidence)
    if args.baseline:
        collectors["baseline"] = baseline_module(args.baseline, "read_collector").ReadCollector()
        extractors["baseline"] = baseline_module(args.baseline, "sv_rna_orfs").record_evidence
    report = dict(version=__version__, python=platform.python_version(), platform=platform.platform(),
                  benchmark_sha256=sha256(Path(__file__).read_bytes()).hexdigest(),
                  bundle_sha256=sha256(BUNDLE.read_bytes()).hexdigest(), groups={})
    if args.baseline:
        report["baseline_read_collector_sha256"] = sha256(
            (args.baseline / "isovar/read_collector.py").read_bytes()).hexdigest()
    for technology, queries in groups.items():
        rows = {}
        for compact in (False, True):
            def collect(collector):
                return [collector.collect_locus_reads_with_optional_compaction(
                    alignment_file=reads, compact=compact, **kwargs) for reads, kwargs in queries]
            results, timing = measure({name: lambda c=c: collect(c) for name, c in collectors.items()}, args.repeats)
            for name, result in results.items():
                # Fingerprint includes every public/compact observation field.
                serial = [[{field: getattr(read, field) for field in read.__slots__} for read in query]
                          for query in result]
                timing[name]["evidence_sha256"] = sha256(json.dumps(serial, sort_keys=True, default=list).encode()).hexdigest()
                timing[name]["retained_reads"] = sum(map(len, result))
            rows["compact" if compact else "public"] = timing
        reads = [r for records, _ in queries for r in records.reads]
        _, rows["metadata"] = measure({name: lambda f=f: [f(r) for r in reads]
                                       for name, f in extractors.items()}, args.repeats)
        report["groups"][technology] = dict(queries=len(queries), records=len(reads),
            query_bases=sum(r.query_length or 0 for r in reads), measurements=rows)
    with args.output.open("x") as handle:
        json.dump(report, handle, indent=2)
        handle.write("\n")


if __name__ == "__main__":
    main()
