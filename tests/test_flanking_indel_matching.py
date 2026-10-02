"""Anchored overlap costs and spans retain indels without inventing coverage."""

from itertools import product
from types import SimpleNamespace

import pytest

from isovar import AlleleRead
from isovar.rna_candidate_support import _comparison, _flank_alignment

TX = SimpleNamespace(id="tx", exon_intervals=[(1, 1000)])


def compare(read_flank, candidate_flank, side="prefix", allele="G", alt="G", ref="A", **kwargs):
    # Inputs run outward from the focal allele on either side.
    prefix, suffix = (read_flank[::-1], "TTTAAA") if side == "prefix" else ("AAACCC", read_flank)
    before, after = ((candidate_flank[::-1], "TTTAAA") if side == "prefix"
                     else ("AAACCC", candidate_flank))
    read = AlleleRead(prefix, allele, suffix, "read")
    candidate = dict(candidate_id="a", sequence=before + alt + after,
                     variant_interval=[len(before), len(before) + len(alt)])
    return _comparison(read, candidate, ref, [TX], kwargs.get("max_edits", 1), kwargs.get("min_overlap", 1))


@pytest.mark.parametrize("side", ["prefix", "suffix"])
def test_one_internal_insertion_costs_one_edit_not_two(side):
    result = compare("CATTGCA", "CATGCA", side)
    assert result["edit_distance"] == 1
    assert result["reference_edit_distance"] == 2
    assert result["within_tolerance"] and not result["reference_preferred"]
    assert result["candidate_interval"] == [0, 13]
    assert result["spans_candidate"]
    assert not compare("CATTGCA", "CATGCA", side, max_edits=0)["within_tolerance"]


@pytest.mark.parametrize("side", ["prefix", "suffix"])
@pytest.mark.parametrize("position", [0, 1, 5, 11])
@pytest.mark.parametrize("operation", ["insertion", "deletion"])
def test_anchor_and_internal_indels_are_charged_but_distal_overhang_is_free(side, position, operation):
    sequence = "ACGTACGTTAGC"
    read = (sequence[:position] + "G" + sequence[position:] if operation == "insertion"
            else sequence[:position] + sequence[position + 1:])
    expected = 0 if operation == "deletion" and position == len(sequence) - 1 else 1
    assert compare(read, sequence, side)["edit_distance"] == expected
    # Swapping read/candidate retains the same anchored overlap cost.
    assert compare(sequence, read, side)["edit_distance"] == expected


@pytest.mark.parametrize("side", ["prefix", "suffix"])
def test_partial_windows_and_outer_extensions_do_not_add_boundary_errors(side):
    assert compare("ACGT", "ACGTACGTTAGC", side)["edit_distance"] == 0
    assert compare("ACGTACGTTAGC", "ACGT", side)["edit_distance"] == 0
    assert compare("ACTGT", "ACGTACGTTAGC", side)["edit_distance"] == 1
    assert compare("ACGACGTTAGC", "ACGTACGTTAGC", side)["edit_distance"] == 1
    assert compare("ACGTNNN", "ACGT", side)["edit_distance"] == 0
    assert compare("ACGT", "ACGTNNN", side)["edit_distance"] == 0
    assert compare("ACNT", "ACGT", side) is None


@pytest.mark.parametrize("side", ["prefix", "suffix"])
@pytest.mark.parametrize("allele,alt,ref", [("A", "G", "A"), ("T", "G", "A"),
                                           ("", "GG", ""), ("AA", "", "AA")])
def test_flanking_indel_cannot_rescue_reference_or_tied_focal_evidence(side, allele, alt, ref):
    result = compare("CATTGCA", "CATGCA", side, allele=allele, alt=alt, ref=ref, max_edits=10)
    assert result["reference_preferred"]


@pytest.mark.parametrize("alt,ref", [("GG", ""), ("", "AA")])
def test_focal_indels_remain_globally_anchored(alt, ref):
    result = compare("CATTGCA", "CATGCA", allele=alt, alt=alt, ref=ref)
    assert result["edit_distance"] == 1
    assert result["reference_edit_distance"] == 3
    assert not result["reference_preferred"]


def test_overlap_alignment_against_independent_dynamic_programming():
    # Enumerate a full unit-cost DP matrix. Valid stopping points are the last
    # row/column, with neither start free. This does not call Edlib.
    words = ["".join(letters) for n in range(5) for letters in product("AC", repeat=n)]
    for read in words:
        for candidate in words:
            n, m = len(read), len(candidate)
            costs = [[i + j for j in range(m + 1)] for i in range(n + 1)]
            for i in range(1, n + 1):
                for j in range(1, m + 1):
                    costs[i][j] = min(costs[i - 1][j] + 1, costs[i][j - 1] + 1,
                                      costs[i - 1][j - 1] + (read[i - 1] != candidate[j - 1]))
            boundary = {(n, j) for j in range(m + 1)} | {(i, m) for i in range(n + 1)}
            best = min(costs[i][j] for i, j in boundary)
            endpoints = [list(pair) for pair in sorted(boundary) if costs[pair[0]][pair[1]] == best]
            assert _flank_alignment(read, candidate) == dict(edit_distance=best, endpoints=endpoints)


def test_equal_lengths_still_require_both_alignment_directions():
    # A read insertion followed by a candidate overhang is cheaper than
    # consuming the full candidate against the read's prefix.
    assert _flank_alignment("AAC", "CAA") == dict(edit_distance=1, endpoints=[[2, 3]])
    assert _flank_alignment("CAA", "AAC") == dict(edit_distance=1, endpoints=[[3, 2]])


@pytest.mark.parametrize("side", ["prefix", "suffix"])
def test_tied_endpoints_do_not_assert_full_coverage_or_meet_a_larger_overlap_floor(side):
    result = compare("ACGT", "ACGA", side, min_overlap=10)
    assert result["edit_distance"] == 1
    assert result["flank_alignments"][side]["endpoints"] == [[3, 4], [4, 3], [4, 4]]
    assert result["compared_read_bases"] == 10  # 3 certain flank bases + allele + 6 exact other-side bases
    assert not result["spans_candidate"]
    assert result["candidate_interval"] == ([1, 11] if side == "prefix" else [0, 10])
    assert compare("ACGT", "ACGA", side, min_overlap=11) is None


def test_empty_flanks_add_neither_coverage_nor_path_requirements():
    assert _flank_alignment("", "NNN") == dict(edit_distance=0, endpoints=[[0, 0]])
    assert _flank_alignment("NNN", "") == dict(edit_distance=0, endpoints=[[0, 0]])
    result = compare("", "NNN", min_overlap=7)
    assert result["edit_distance"] == 0 and result["compared_read_bases"] == 7
    assert not result["spans_candidate"]
    assert compare("", "NNN", min_overlap=8) is None


def test_path_checks_include_the_longest_tied_read_span():
    candidate = dict(candidate_id="a", sequence="GACGA", variant_interval=[0, 1])
    transcript = SimpleNamespace(id="tx", exon_intervals=[(1, 4), (101, 110)])
    # The final base, included by one optimal endpoint, crosses to a different
    # exon without a recorded splice. Taking the shortest endpoint would hide it.
    read = AlleleRead("", "G", "ACGT", "read", reference_blocks=((0, 4, 0, 4), (4, 5, 100, 101)))
    assert _comparison(read, candidate, "A", [transcript], 1, 4) is None
    spliced = AlleleRead("", "G", "ACGT", "read", reference_blocks=read.reference_blocks,
                         splice_junctions=((4, 100),))
    assert _comparison(spliced, candidate, "A", [transcript], 1, 4)["edit_distance"] == 1
    unknown = AlleleRead("", "G", "ACGN", "read")
    assert _comparison(unknown, candidate, "A", [transcript], 1, 4) is None


def test_unobserved_distal_path_is_neutral_after_an_insertion():
    candidate = dict(candidate_id="a", sequence="GCATGCA", variant_interval=[0, 1])
    # Candidate end is reached after 7 query flank bases (the inserted T counts),
    # so the unrelated path of the last two read bases is outside the comparison.
    read = AlleleRead("", "G", "CATTGCAGG", "read",
                      reference_blocks=((0, 8, 0, 8), (8, 10, 100, 102)), splice_junctions=((8, 100),))
    result = _comparison(read, candidate, "A", [TX], 1, 8)
    assert result["edit_distance"] == 1 and result["compared_read_bases"] == 8
    assert result["spans_candidate"]
