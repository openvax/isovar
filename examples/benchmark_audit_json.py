"""Measure an audit JSON stage in a fresh process against the released helpers.

Use --trace separately for incremental allocations; tracing affects timings.
Parsing is measured with process RSS, since tracking every decoded object would
itself add substantial memory. Input bytes are always checksum-pinned.
"""

import argparse
import gc
import gzip
from hashlib import sha256
import json
from pathlib import Path
import resource
import sys
import tempfile
import time
import tracemalloc

from isovar._audit_json import canonical, identity, normalize_json
from isovar.audit_bundle import digest
from examples.sid_sv_audit.inventory import read_json, write_json


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("input", type=Path)
    parser.add_argument("--operation", choices=("read", "hash", "write", "normalize"), required=True)
    parser.add_argument("--implementation", choices=("old", "new"), required=True)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--trace", action="store_true")
    args = parser.parse_args()
    if args.trace and args.operation == "read":
        parser.error("Measure parsing separately without allocation tracing")
    source_sha256 = digest(args.input)
    old = args.implementation == "old"
    value = None if args.operation == "read" else read_json(args.input)
    gc.collect()
    with tempfile.TemporaryDirectory(dir=args.output.parent) as directory:
        path = Path(directory) / "result.json.gz"
        if args.trace:
            tracemalloc.start()
        start, cpu = time.perf_counter(), time.process_time()
        if args.operation == "read":
            if old:
                data = args.input.read_bytes()
                value = json.loads(gzip.decompress(data) if args.input.suffix == ".gz" else data)
                del data
            else:
                value = read_json(args.input)
        elif args.operation == "hash":
            checksum = sha256(canonical(value)).hexdigest() if old else identity(value)
        elif args.operation == "write":
            if old:
                path.write_bytes(gzip.compress(canonical(value) + b"\n", mtime=0))
            else:
                write_json(path, value)
        else:
            normalized = json.loads(canonical(value)) if old else normalize_json(value)
        elapsed, cpu_elapsed = time.perf_counter() - start, time.process_time() - cpu
        peak_allocation = tracemalloc.get_traced_memory()[1] if args.trace else None
        if args.trace:
            tracemalloc.stop()
        rss = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss * (1 if sys.platform == "darwin" else 1024)
        report = dict(operation=args.operation, implementation=args.implementation,
                      input_sha256=source_sha256, input_bytes=args.input.stat().st_size,
                      seconds=elapsed, cpu_seconds=cpu_elapsed, peak_rss_bytes=rss,
                      traced_peak_incremental_bytes=peak_allocation)
        if args.operation == "hash":
            report["canonical_sha256"] = checksum
        elif args.operation == "write":
            report.update(output_sha256=digest(path), output_bytes=path.stat().st_size)
        else:
            report["canonical_sha256"] = identity(normalized if args.operation == "normalize" else value)
    args.output.write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps(report), flush=True)


if __name__ == "__main__":
    main()
