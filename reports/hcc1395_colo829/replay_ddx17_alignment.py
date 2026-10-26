"""Replay STAR on the original DDX17 regional pairs and selected-pair rescue."""

import argparse
import gzip
from pathlib import Path
import subprocess

from .build_ddx17_audit import HERE


def replay(star, work_dir):
    if subprocess.check_output([star, '--version'], text=True).strip() != '2.7.11b':
        raise ValueError('This replay pins STAR 2.7.11b')
    work_dir.mkdir(parents=True, exist_ok=True)
    index = work_dir / 'index'
    index.mkdir(exist_ok=True)
    pairs = {}
    for name in ('regional-pairs.fastq.gz', 'selected-pairs.fastq.gz'):
        with gzip.open(HERE / name, 'rt') as stream:
            while header := stream.readline():
                fragment, mate = header.split()[0].rsplit('.', 1)
                record = [fragment + '\n', stream.readline(), stream.readline(), stream.readline()]
                if not record[-1]:
                    raise ValueError('Truncated FASTQ')
                key = fragment, mate
                if key in pairs and pairs[key] != record:
                    raise ValueError(f'Conflicting original read: {key}')
                pairs[key] = record
    if len(pairs) != 27309 * 2:
        raise ValueError('Unexpected pair inventory')
    fastqs = [work_dir / f'reads.{mate}.fastq' for mate in ('1', '2')]
    for mate, path in zip(('1', '2'), fastqs):
        with path.open('w') as stream:
            for key in sorted(pairs):
                if key[1] == mate:
                    stream.writelines(pairs[key])
    subprocess.run([star, '--runMode', 'genomeGenerate', '--genomeDir', str(index),
                    '--genomeFastaFiles', str(HERE / 'reference.fa'),
                    '--genomeSAindexNbases', '6', '--genomeChrBinNbits', '14',
                    '--runThreadN', '2', '--outFileNamePrefix', str(work_dir / 'index.')], check=True)
    subprocess.run([star, '--genomeDir', str(index), '--readFilesIn', *map(str, fastqs),
                    '--runThreadN', '2', '--outFileNamePrefix', str(work_dir / 'expanded.star.'),
                    '--outSAMtype', 'BAM', 'SortedByCoordinate', '--outReadsUnmapped', 'Fastx',
                    '--outSAMattributes', 'Standard', '--alignIntronMax', '50000'], check=True)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--star', required=True)
    parser.add_argument('--work-dir', type=Path, required=True)
    args = parser.parse_args()
    replay(args.star, args.work_dir.resolve())
