"""Sequence continuity and original-spot accounting for SV RNA paths."""

from isovar.sv_rna_path_support import (
    PathMatch, exact_path_matches, interval_spots, overlap_chain, reverse_complement)


def test_exact_path_support_deduplicates_alternative_mappings_and_tracks_clips():
    path = "AGTTCGACCTGATCGTACGATTCAGGCTACCGTAAGTCTGACCTATGCGTTAACGATC"
    forward = path[5:45]
    reverse = reverse_complement(path[25:60])
    reads = [
        ("spot-a", 1, forward, 0, 25),
        ("spot-a", 1, forward, 0, 25),  # Another alignment of the same end.
        ("spot-b", 2, reverse, None, None),
    ]
    matches = exact_path_matches(path, reads)
    assert [match for match in matches if match.spot == "spot-a"] == [
        PathMatch("spot-a", 1, 5, 45, clipped=True)]
    assert PathMatch("spot-b", 2, 25, len(path)) in matches
    assert interval_spots(matches, 30, 38) == ["spot-a", "spot-b"]


def test_mate_pair_is_not_an_exact_sequence_bridge():
    first = PathMatch("one-spot", 1, 0, 30)
    mate = PathMatch("one-spot", 2, 20, 50)
    assert overlap_chain([first, mate], 5, 45, min_overlap=10) == []
    independent = PathMatch("second-spot", 2, 20, 50)
    assert overlap_chain([first, independent], 5, 45, min_overlap=10) == [
        dict(spot="one-spot", mate=1, path_interval=[0, 30]),
        dict(spot="second-spot", mate=2, path_interval=[20, 50]),
    ]
    assert overlap_chain([first, independent], 5, 45, min_overlap=11) == []


def test_overlap_chain_revisits_interval_with_different_fragment_history():
    from isovar.sv_rna_path_support import overlap_chain_search
    matches = [PathMatch('a', 1, 0, 10, False), PathMatch('b', 1, 0, 9, False),
               PathMatch('c', 1, 5, 15, False), PathMatch('a', 2, 10, 20, False)]
    result = overlap_chain_search(matches, 0, 20, min_overlap=4, max_reads=3)
    assert result['status'] == 'found'
    assert [m['spot'] for m in result['chain']] == ['b', 'c', 'a']
    limited = overlap_chain_search(matches, 0, 20, min_overlap=4, max_reads=3, max_states=1)
    assert limited['status'] == 'search_limit' and not limited['chain']


def test_dense_equivalent_intervals_are_bounded_and_gaps_reject_without_search():
    from isovar.sv_rna_path_support import overlap_chain_search
    matches = [PathMatch(str(i), 1, a, a + 100, False) for a in range(0, 901, 50) for i in range(1000)]
    result = overlap_chain_search(matches, 0, 1000, min_overlap=30, max_reads=20)
    assert result['status'] == 'found' and len(result['chain']) == 19
    assert len({m['spot'] for m in result['chain']}) == 19
    assert result['states_examined'] == 19
    gap = overlap_chain_search(matches, 0, 1100, min_overlap=30, max_reads=20)
    assert gap['status'] == 'no_chain_within_read_limit' and gap['states_examined'] == 0


def test_chain_search_agrees_with_independent_small_exhaustive_enumeration():
    from itertools import permutations
    import random
    from isovar.sv_rna_path_support import overlap_chain_search
    randomizer = random.Random(13)
    for _ in range(100):
        matches = [PathMatch(str(randomizer.randrange(4)), 1, a := randomizer.randrange(9),
                             a + randomizer.randrange(1, 6), False) for _ in range(8)]
        minimum, depth = randomizer.randrange(1, 4), randomizer.randrange(1, 5)
        found = False
        for length in range(1, depth + 1):
            for chain in permutations(matches, length):
                if len({m.spot for m in chain}) != length or chain[0].start > 0:
                    continue
                reach = chain[0].end
                for match in chain[1:]:
                    if match.start > reach - minimum or match.end <= reach:
                        break
                    reach = match.end
                else:
                    if reach >= 10:
                        found = True
                        break
            if found:
                break
        result = overlap_chain_search(matches, 0, 10, min_overlap=minimum, max_reads=depth)
        assert found == (result['status'] == 'found')
