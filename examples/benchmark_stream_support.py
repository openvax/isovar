"""Measure support working memory separately from retained synthetic output.

Run identical invocations against each checkout in fresh processes. Tracemalloc
measures Python allocations during reconstruction, not native allocations or
serialization. Peak minus post-GC retained memory is a transient-allocation
estimate, not an independent measurement of per-base geometry alone. RSS includes
native memory. These artificial, error-free reads are not an ONT performance claim.
"""

import argparse
import gc
from hashlib import sha256
import json
from pathlib import Path
import random
import resource
import signal
import sys
import tempfile
import time
import tracemalloc


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--checkout', type=Path, required=True)
    parser.add_argument('--reads', type=int, required=True)
    parser.add_argument('--seconds', type=int, default=600)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    if args.reads < 1 or args.seconds < 1:
        parser.error('Read count and runtime budget must be positive')
    sys.path.insert(0, str(args.checkout.resolve()))
    import pysam
    from isovar.sv_rna import reconstruct_sv_rna
    from isovar.fusion import FusionBreakpoint

    rng = random.Random(42)
    dna = ''.join(rng.choices('ACGT', k=20000))
    header = pysam.AlignmentHeader.from_dict(dict(HD=dict(VN='1.6', SO='coordinate'), SQ=[dict(SN='1', LN=20000)]))

    def expired(*unused):
        raise TimeoutError('Synthetic benchmark runtime budget reached')

    signal.signal(signal.SIGALRM, expired)
    signal.alarm(args.seconds)
    with tempfile.TemporaryDirectory() as directory:
        original, sorted_bam = Path(directory) / 'original.bam', Path(directory) / 'sorted.bam'
        with pysam.AlignmentFile(original, 'wb', header=header) as output:
            for i in reversed(range(args.reads)):
                flank = 1000 + i % 1000
                read = pysam.AlignedSegment(header)
                read.query_name = 'read-%08d' % i
                read.reference_id, read.reference_start, read.mapping_quality = 0, 5000 - flank, 60
                read.cigartuples = [(0, flank), (3, 5000), (0, flank)]
                read.query_sequence = dna[5000 - flank:5000] + dna[10000:10000 + flank]
                read.query_qualities = [30] * (2 * flank)
                output.write(read)
        # Keep the independently generated inputs byte-equivalent: a sort PG
        # command would otherwise embed each process's temporary directory.
        pysam.sort('--no-PG', '-o', str(sorted_bam), str(original))
        pysam.index(str(sorted_bam))
        tracemalloc.start()
        start = time.perf_counter()
        with pysam.AlignmentFile(sorted_bam) as bam:
            result = reconstruct_sv_rna(
                bam, donor=FusionBreakpoint('1', 5000, '+'), acceptor=FusionBreakpoint('1', 10000, '+'),
                event_id='synthetic', reference_name='synthetic', sample_id='sample', source='synthetic',
                event_provenance={'synthetic': True}, references=[], regions=[], dense_support=True,
                max_records=1, max_paths=1, assemble=False, max_orf_candidates=1)
        gc.collect()
        retained, peak = tracemalloc.get_traced_memory()
        tracemalloc.stop()
        elapsed = time.perf_counter() - start
        rss = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss * (1 if sys.platform == 'darwin' else 1024)
        serialized = json.dumps(result, sort_keys=True, separators=(',', ':')).encode()
        report = dict(
            reads=args.reads, seconds_budget=args.seconds, seconds=elapsed,
            driver_sha256=sha256(Path(__file__).read_bytes()).hexdigest(),
            engine_sha256={name: sha256(path.read_bytes()).hexdigest()
                           for name in ('sv_rna.py', 'sv_rna_dense.py', 'sv_rna_stream.py')
                           if (path := args.checkout / 'isovar' / name).exists()},
            python_peak_bytes=peak, retained_python_bytes=retained, transient_python_bytes=peak - retained,
            peak_rss_bytes=rss, serialized_output_bytes=len(serialized),
            result_sha256=sha256(serialized).hexdigest(), paths=len(result['paths']),
            support=[p['full_interval_support']['reads'] for p in result['paths']])
        args.output.write_text(json.dumps(report, indent=2) + '\n')
    signal.alarm(0)
    print(json.dumps(report))


if __name__ == '__main__':
    main()
