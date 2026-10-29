"""Measure repeated annotation of four unchanged Sid records, not whole-BAM throughput."""

import argparse
from dataclasses import asdict
from pathlib import Path
from importlib.metadata import version
from tempfile import TemporaryDirectory
from time import perf_counter

import pysam

from isovar import __version__, Adapter, ReadEndProfile, ReadEndInferenceCache, sid_data
from isovar._file_identity import file_sha256
from isovar.read_end_inference import read_sequence_view_from_alignment
from isovar.sv_rna_audit import write_json


def benchmark(iterations: int) -> dict:
    if iterations < 1:
        raise ValueError('Iterations must be positive')
    profile = ReadEndProfile('documented-Illumina-R2-motif', adapters=(Adapter('R2',
        'AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT', mate=2,
        source='https://github.com/Illumina/isaac2/blob/master/src/markdown/manual.md'),))
    records, inputs = [], []
    with TemporaryDirectory(prefix='isovar-read-end-benchmark-') as folder:
        for index, member in enumerate(sorted(name for name in sid_data.members()
                                              if name.startswith('read_ends/osteosarc.sam#'))):
            path = sid_data.export(member, Path(folder) / f'{index}.bam')
            with pysam.AlignmentFile(path) as bam:
                read, = list(bam)
                records.append(read)
                inputs.append(dict(member=member, exported_bam_sha256=file_sha256(path),
                                   name=read.query_name, flag=read.flag, cigar=read.cigarstring))
    originals = [read.to_string() for read in records]
    baseline = tuple(read_sequence_view_from_alignment(read, profile, trim_adapters=True, trim_poly_a=True)
                     for read in records)
    timings = {}
    cache = ReadEndInferenceCache()
    for name, annotation_cache in (('uncached', None), ('cached', cache)):
        started = perf_counter()
        for _ in range(iterations):
            result = tuple(read_sequence_view_from_alignment(read, profile,
                trim_adapters=True, trim_poly_a=True, cache=annotation_cache) for read in records)
            if result != baseline:
                raise AssertionError('Cached annotation changed evidence or sequence eligibility')
        timings[name] = perf_counter() - started
    assert originals == [read.to_string() for read in records]
    return dict(schema='isovar.read_end_cache_benchmark.v1', producer_version=__version__,
        implementation_sha256={name: file_sha256(Path(__file__).resolve().parents[1] / 'isovar' / name)
                               for name in ('read_end_inference.py', 'read_end_profiles.py', 'read_collector.py')},
        dependencies={name: version(name) for name in ('pysam', 'edlib')}, driver_sha256=file_sha256(__file__),
        scope='Four unchanged selected original Sid records; repeated annotation only, not whole-library throughput',
        inputs=inputs, profile=asdict(profile), profile_sha256=profile.fingerprint,
        iterations=iterations, view_requests=iterations * len(records), seconds=timings,
        cache=dict(hits=cache.hits, misses=cache.misses, retained_bases=cache.retained_bases,
                   max_bases=cache.max_bases, max_entries=cache.max_entries),
        unchanged_originals=True, identical_derived_views=True,
        limitation='A documented R2 motif does not identify the actual Sid library preparation; no false fusion is claimed')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--iterations', type=int, default=200)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    args.output.parent.mkdir(parents=True, exist_ok=True)
    write_json(args.output, benchmark(args.iterations))
