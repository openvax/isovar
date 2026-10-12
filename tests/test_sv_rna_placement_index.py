import random
from types import SimpleNamespace

import pytest

from isovar.sv_rna import _merge
from isovar.sv_rna_dense import _overlaps_placements
from isovar.sv_rna_orfs import _PlacementIndex, _witnesses


def test_run_placement_comparison_matches_literal_per_base_rule():
    rng = random.Random(470)

    def positions():
        values = []
        for _ in range(rng.randrange(1, 12)):
            low, size, strand, contig = rng.randrange(50), rng.randrange(1, 15), rng.choice('+-'), rng.choice('12')
            bases = range(low, low + size) if strand == '+' else range(low + size - 1, low - 1, -1)
            values.extend((contig, q, strand) for q in bases)
            values.extend([None] * rng.randrange(5))
        return tuple(values)

    for case in range(250):
        own = positions()
        other = (None,) * 4 + own + (None,) * 3 if case % 2 else positions()
        left, right = _PlacementIndex(own), _PlacementIndex(other)
        for shift in (-7, -3, 0, 4, 8):
            lo, hi = max(0, -shift), min(len(own), len(other) - shift)
            if hi < lo:
                continue
            start = rng.randrange(lo, hi + 1)
            end = rng.randrange(start, hi + 1)
            expected = all(p is None or q is None or p == q
                           for p, q in zip(own[start:end], other[start + shift:end + shift]))
            assert left.agrees(right, start, end, shift) == expected


@pytest.mark.parametrize('blocks,expected', [
    ([(0, 10)], False), ([(30, 40)], False), ([(30, 41)], True),
    ([(11, 20)], True), ([(0, 10), (60, 70)], False), ([], False),
])
def test_candidate_union_prefilter_keeps_cigar_gaps_and_half_open_boundaries(blocks, expected):
    intervals = _merge([(10, 20), (20, 30), (40, 50), (45, 60)])
    read = SimpleNamespace(get_blocks=lambda: blocks)
    assert bool(_overlaps_placements(read, intervals)) == expected
    assert not _overlaps_placements(read, [])


def test_witness_cache_exhaustion_keeps_repeated_occurrences_and_links(monkeypatch):
    from isovar import sv_rna_orfs

    monkeypatch.setattr(sv_rna_orfs, '_WITNESS_CACHE_ENTRIES', 1)
    positions = tuple(('1', q, '+') for q in range(6))
    observed = (None,) * 2 + positions + (None,) * 2 + positions
    observation = SimpleNamespace(key='original', positions=observed, sequence='A' * 16,
                                  query_interval=(0, 16), reverse_complement=False,
                                  comparison_key=('A' * 16, observed))
    junction = dict(query_interval=[2, 2], left=list(positions[1]), right=list(positions[2]),
                    direct_observations=['original'])
    witnesses = _witnesses('A' * 6, positions, 1, 3, [junction, junction], {'original': observation})
    assert [row['original_query_interval'] for row in witnesses] == [[3, 5], [11, 13]]
    assert all(len(row['junction_links']) == 1 for row in witnesses)


@pytest.mark.parametrize('entries,bases', [(4096, 250_000), (1, 250_000), (4096, 1)])
def test_shared_exact_search_keeps_junction_membership_and_repeated_witnesses(monkeypatch, entries, bases):
    from isovar import sv_rna_orfs

    monkeypatch.setattr(sv_rna_orfs, '_WITNESS_CACHE_ENTRIES', entries)
    monkeypatch.setattr(sv_rna_orfs, '_WITNESS_CACHE_BASES', bases)
    positions = tuple(('1', q, '+') for q in range(6))
    observed = (None,) * 2 + positions + (None,) * 2 + positions
    observations = {
        name: SimpleNamespace(key=name, positions=observed, sequence='A' * 16,
                              query_interval=(0, 16), reverse_complement=False,
                              comparison_key=('A' * 16, observed))
        for name in ('first', 'second', 'unrelated')}
    junction = dict(query_interval=[2, 2], left=list(positions[1]), right=list(positions[2]),
                    direct_observations=['first'])
    other = dict(junction, direct_observations=['second'])
    search = sv_rna_orfs._SequenceMatches(o.sequence for o in observations.values())
    cache = {sv_rna_orfs._SEQUENCE_MATCHES: search}
    for start, end in ((1, 3), (2, 4)):
        got = _witnesses('A' * 6, positions, start, end, [junction, other], observations, group_cache=cache)
        expected = _witnesses('A' * 6, positions, start, end, [junction, other], observations)
        assert got == expected
        assert [row['observation'] for row in got] == ['first', 'first', 'second', 'second']
        assert all(len(row['junction_links']) == 1 for row in got)
    assert search.entries <= entries and search.bases <= bases
