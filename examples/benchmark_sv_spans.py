"""Compare exact candidate spans against a preceding checkout's implementation.

This isolates fragmented-path matching; it does not measure reconstruction.
Run from the current checkout, supplying the preceding checkout as --baseline.
"""

import argparse
import ast
import json
from pathlib import Path
import random
import time
from types import SimpleNamespace

from isovar import sv_rna


def compare(baseline, comparisons=1000):
    source = (baseline / 'isovar/sv_rna.py').read_text()
    node = next(n for n in ast.parse(source).body if isinstance(n, ast.ClassDef) and n.name == '_PathSpans')
    namespace = dict(vars(sv_rna))
    exec(compile(ast.Module(body=[node], type_ignores=[]), str(baseline), 'exec'), namespace)
    previous = namespace['_PathSpans']
    rng = random.Random(485)

    def positions():
        values = []
        for _ in range(rng.randint(1, 15)):
            contig, start, strand, size = rng.choice(['1', '2']), rng.randrange(200), rng.choice('+-'), rng.randrange(1, 30)
            bases = range(start, start + size) if strand == '+' else range(start + size - 1, start - 1, -1)
            values.extend((contig, base, strand) for base in bases)
            values.extend([None] * rng.randrange(4))
        return sv_rna._Positions(values)

    for index in range(1000):
        pos = positions()
        sequence = ''.join(rng.choices('ACGT', k=len(pos)))
        other = pos if index % 3 == 0 else positions()
        observation = SimpleNamespace(key=str(index), positions=other,
                                      sequence=sequence if other == pos else ''.join(rng.choices('ACGT', k=len(other))))
        assert previous(sequence, pos)(observation) == sv_rna._PathSpans(sequence, pos)(observation), index

    pos = sv_rna._Positions([('1', q + q // 10, '+') for q in range(4000)])
    sequence = 'ACGT' * 1000
    observation = SimpleNamespace(positions=pos, sequence=sequence)
    times = []
    for implementation in (previous, sv_rna._PathSpans):
        matcher = implementation(sequence, pos)
        start = time.perf_counter()
        for index in range(comparisons):
            observation.key = str(index)
            assert matcher(observation) == (0, 3999)
        times.append(time.perf_counter() - start)
    return dict(random_cases=1000, random_equal=True, fragmented_comparisons=comparisons,
                reference_seconds=times[0], indexed_seconds=times[1], component_speedup=times[0] / times[1])


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--baseline', type=Path, required=True)
    parser.add_argument('--comparisons', type=int, default=1000)
    args = parser.parse_args()
    if args.comparisons < 1:
        parser.error('--comparisons must be positive')
    print(json.dumps(compare(args.baseline, args.comparisons), indent=2))
