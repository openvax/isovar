"""Measure one pinned, local audit-bundle view with a selected engine checkout.

Run this identical driver in fresh processes against each checkout and budget.
Times cover reconstruction, including spooling, grouping and complete support;
input verification is outside the timed interval. No remote input is acquired.
"""

import argparse
import cProfile
import gzip
from hashlib import sha256
import json
import os
from pathlib import Path
import platform
import resource
import signal
import socket
import sys
import tempfile
import threading
import time


def checksum(path):
    value = sha256()
    with path.open('rb') as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b''):
            value.update(block)
    return value.hexdigest()


def fingerprint(value):
    return sha256(json.dumps(value, sort_keys=True, separators=(',', ':')).encode()).hexdigest()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--checkout', required=True, type=Path)
    parser.add_argument('--bundle', required=True, type=Path)
    parser.add_argument('--run', required=True)
    parser.add_argument('--max-records', required=True, type=int)
    parser.add_argument('--seconds', type=int, default=600)
    parser.add_argument('--scratch-dir', type=Path)
    parser.add_argument('--output', required=True, type=Path)
    parser.add_argument('--profile', type=Path, help='Profile reconstruction only; excludes report serialization')
    parser.add_argument('--max-scratch-bytes', type=int, help='Abort if sampled temporary storage exceeds this budget')
    args = parser.parse_args()
    if args.max_records < 1 or args.seconds < 1 or (args.max_scratch_bytes is not None and args.max_scratch_bytes < 1):
        parser.error('Budgets must be positive')
    sys.path.insert(0, str(args.checkout.resolve()))
    import pysam
    import isovar
    from isovar import sv_rna_dense
    from isovar.audit_bundle import AuditBundle
    from isovar.audit_replay import read_json
    from isovar.read_collector import ReadCollector
    from isovar.sv_rna import reconstruct_sv_rna, sv_rna_input_from_dict

    bundle = AuditBundle.open(args.bundle, verify=False)
    plan = bundle.manifest['metadata']['sv_rna_runs'][args.run]
    pins = {}
    for key in ('bam', 'index', 'input'):
        path = bundle.path(plan[key])
        expected = bundle.manifest['files'][plan[key]]
        assert path.stat().st_size == expected['size'] and checksum(path) == expected['sha256']
        pins[key] = expected
    inputs = sv_rna_input_from_dict(read_json(bundle.path(plan['input'])))
    parameters = dict(plan['parameters'], dense_support=True, max_records=args.max_records)
    collector = ReadCollector(**plan.get('read_collection', {}))
    stages, scratch_peak, scratch_exceeded, stopped = {}, [0], [False], threading.Event()

    def timed(name, function):
        def run(*values, **kwargs):
            start = time.perf_counter()
            try:
                return function(*values, **kwargs)
            finally:
                stages[name] = time.perf_counter() - start
        return run

    for name in ('_spool', '_select', '_support'):
        setattr(sv_rna_dense, name, timed(name, getattr(sv_rna_dense, name)))

    def offline(*args, **kwargs):
        raise AssertionError('Benchmark attempted network access')

    socket.socket.connect = offline
    socket.create_connection = offline

    def timeout(*unused):
        raise TimeoutError('Declared reconstruction time budget reached')

    signal.signal(signal.SIGALRM, timeout)
    report = dict(version=isovar.__version__, python=platform.python_version(), platform=platform.platform(),
                  driver_sha256=checksum(Path(__file__)), bundle_id=bundle.bundle_id, run_id=args.run,
                  input_files=pins, parameters=parameters, read_collection=plan.get('read_collection', {}),
                  seconds_budget=args.seconds, scratch_bytes_budget=args.max_scratch_bytes,
                  profiled=bool(args.profile), original_input_limitations=plan.get('input_limitations', []),
                  engine_sha256={name: checksum(args.checkout / 'isovar' / name)
                                 for name in ('sv_rna.py', 'sv_rna_dense.py', 'sv_rna_orfs.py',
                                              'read_identity.py', 'read_lineage.py', 'cell_umi.py')})
    with tempfile.TemporaryDirectory(prefix='dense-benchmark-', dir=args.scratch_dir) as scratch:
        def monitor():
            while not stopped.wait(0.2):
                total = 0
                for path in Path(scratch).rglob('*'):
                    try:
                        if path.is_file():
                            total += path.stat().st_size
                    except FileNotFoundError:
                        pass
                scratch_peak[0] = max(scratch_peak[0], total)
                if args.max_scratch_bytes is not None and total > args.max_scratch_bytes:
                    scratch_exceeded[0] = True
                    os.kill(os.getpid(), signal.SIGALRM)
                    return

        thread = threading.Thread(target=monitor, daemon=True)
        thread.start()
        profile = cProfile.Profile() if args.profile else None
        if profile is not None:
            profile.enable()
        start, cpu = time.perf_counter(), time.process_time()
        signal.alarm(args.seconds)
        try:
            with pysam.AlignmentFile(bundle.path(plan['bam']), index_filename=bundle.path(plan['index'])) as bam:
                result = reconstruct_sv_rna(bam, **inputs, source=plan['source'],
                                            read_collector=collector, scratch_dir=scratch, **parameters)
            report.update(status='completed', discovery=result['discovery'],
                          support_acquisition=result['support_acquisition'], limitations=result['limitations'])
        except TimeoutError:
            report['status'] = 'scratch_limit' if scratch_exceeded[0] else 'timeout'
            result = None
        finally:
            signal.alarm(0)
            if profile is not None:
                profile.disable()
            report.update(seconds=time.perf_counter() - start, cpu_seconds=time.process_time() - cpu,
                          stages_seconds=stages, sampled_peak_scratch_bytes=scratch_peak[0],
                          peak_rss_bytes=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
                          * (1 if sys.platform == 'darwin' else 1024))
            stopped.set()
            thread.join()
            if profile is not None:
                profile.dump_stats(args.profile)
    if result is not None:
        report['paths'] = [dict(
            path_id=p['path_id'], sequence_sha256=fingerprint(p['sequence']),
            complete_support_sha256=fingerprint(dict(
                blocks=p['blocks'], junctions=p['junctions'], full_interval_support=p['full_interval_support'],
                orfs=[dict(sequence=o['amino_acids'], query_interval=o['query_interval'],
                           full_interval_support=o['full_interval_support'])
                      for o in p['exploratory_orfs']['candidates']])),
            full_interval_fragments=p['full_interval_support']['fragments']) for p in result['paths']]
        serialized = time.perf_counter()
        with gzip.open(str(args.output) + '.result.json.gz', 'wt') as handle:
            handle.write(json.dumps(result, separators=(',', ':')))
        report['result_serialization_seconds'] = time.perf_counter() - serialized
    args.output.write_text(json.dumps(report, indent=2) + '\n')
    print(json.dumps({k: report[k] for k in ('status', 'seconds', 'cpu_seconds', 'peak_rss_bytes', 'stages_seconds')}))
    return 0 if report['status'] == 'completed' else 2


if __name__ == '__main__':
    raise SystemExit(main())
