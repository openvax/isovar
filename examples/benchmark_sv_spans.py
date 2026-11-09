"""Compare exact candidate spans against a preceding checkout's implementation.

This isolates fragmented-path matching; it does not measure reconstruction.
Run from the current checkout, supplying the preceding checkout as --baseline.
"""

import argparse
import ast
import json
import random
import time
from pathlib import Path
from types import SimpleNamespace

import isovar.rna_path_geometry as rna_path_geometry


def compare(baseline, comparisons=1000):
    source_path = baseline / 'isovar/rna_path_geometry.py'
    inclusive_baseline = not source_path.exists()
    if inclusive_baseline:
        source_path = baseline / 'isovar/sv_rna.py'
    source = source_path.read_text()
    class_name = '_PathSpans' if inclusive_baseline else 'PathSpans'
    node = next(n for n in ast.parse(source).body if isinstance(n, ast.ClassDef) and n.name == class_name)
    namespace = dict(vars(rna_path_geometry), _runs=rna_path_geometry.position_runs)
    exec(compile(ast.Module(body=[node], type_ignores=[]), str(baseline), 'exec'), namespace)
    previous = namespace[class_name]

    def span(matcher, observation, inclusive):
        result = matcher(observation)
        return (result[0], result[1] + 1) if result and inclusive else result
    rng = random.Random(485)

    def positions():
        values = []
        for _ in range(rng.randint(1, 15)):
            contig, start, strand, size = rng.choice(['1', '2']), rng.randrange(200), rng.choice('+-'), rng.randrange(1, 30)
            bases = range(start, start + size) if strand == '+' else range(start + size - 1, start - 1, -1)
            values.extend((contig, base, strand) for base in bases)
            values.extend([None] * rng.randrange(4))
        return rna_path_geometry.PathPositions(values)

    for index in range(1000):
        pos = positions()
        sequence = ''.join(rng.choices('ACGT', k=len(pos)))
        other = pos if index % 3 == 0 else positions()
        observation = SimpleNamespace(key=str(index), positions=other,
                                      sequence=sequence if other == pos else ''.join(rng.choices('ACGT', k=len(other))))
        assert span(previous(sequence, pos), observation, inclusive_baseline) == rna_path_geometry.PathSpans(sequence, pos)(observation), index

    pos = rna_path_geometry.PathPositions([('1', q + q // 10, '+') for q in range(4000)])
    sequence = 'ACGT' * 1000
    observation = SimpleNamespace(positions=pos, sequence=sequence)
    times = []
    for implementation, inclusive in ((previous, inclusive_baseline), (rna_path_geometry.PathSpans, False)):
        matcher = implementation(sequence, pos)
        start = time.perf_counter()
        for index in range(comparisons):
            observation.key = str(index)
            assert span(matcher, observation, inclusive) == (0, 4000)
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
