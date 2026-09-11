"""Package #229 collection and whole-pipeline comparisons, failing closed."""

import argparse
from pathlib import Path

from tests.data.osteosarc.expansion.inventory import digest, write_json
from tests.data.osteosarc.expansion.performance_report import collect_benchmarks, read_json


def collect_runs(entries, comparisons):
    runs = {}
    for entry in entries:
        label, directory = entry.split("=", 1)
        if label in runs:
            raise ValueError("Duplicate collection benchmark label")
        directory = Path(directory)
        run = dict(identity=read_json(directory / "identity.json"), result=read_json(directory / "result.json"))
        if (directory / "allocations.json").exists():
            run["allocations"] = read_json(directory / "allocations.json")
        run["artifacts"] = {p.name: digest(p) for p in sorted(directory.iterdir()) if p.is_file()}
        runs[label] = run
    checked = []
    for pair in comparisons:
        before, after = pair.split("=", 1)
        left, right = runs[before], runs[after]
        # Pin the measuring script and original acquisition as well as inputs:
        # runs from different benchmark versions are not a valid comparison.
        for key in ("source_bam_sha256", "receipt_sha256", "input_bam_sha256", "input_index_sha256",
                    "inventory_sha256", "benchmark_sha256", "variant_id", "mode",
                    "merge_overlapping_fragments", "allocation_instrumented"):
            if left["identity"][key] != right["identity"][key]:
                raise ValueError(f"Different collection inputs/instrumentation: {pair} {key}")
        if any(run["result"]["status"] != "ok" for run in (left, right)):
            raise ValueError(f"Incomplete collection comparison: {pair}")
        for key in ("locus_reads_sha256", "allele_reads_sha256", "locus_read_count", "allele_read_count", "counts"):
            if left["result"][key] != right["result"][key]:
                raise ValueError(f"Collection results changed: {pair} {key}")
        checked.append(dict(before=before, after=after, exact_locus_and_allele_reads=True,
                            allocation_instrumented=left["identity"]["allocation_instrumented"]))
    return dict(runs=runs, comparisons=checked)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--collection", action="append", default=[], help="label=directory")
    parser.add_argument("--compare-collection", action="append", default=[])
    parser.add_argument("--pipeline", action="append", default=[], help="label=directory")
    parser.add_argument("--compare-pipeline", action="append", default=[])
    parser.add_argument("--output", required=True)
    args = parser.parse_args()
    collection = collect_runs(args.collection, args.compare_collection)
    pipeline = collect_benchmarks(args.pipeline, args.compare_pipeline)
    output = Path(args.output)
    output.mkdir(parents=True, exist_ok=False)
    write_json(output / "collection.json", collection)
    write_json(output / "pipeline.json", pipeline)
    write_json(output / "manifest.json", dict(
        files={name: digest(output / name) for name in ("collection.json", "pipeline.json")},
        generator_sha256=digest(__file__)))
