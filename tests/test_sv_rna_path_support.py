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
