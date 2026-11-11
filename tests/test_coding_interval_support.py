"""Partial contributors and merged mates are not single-read coding spans."""

import pytest

from isovar.allele_read import AlleleRead
from isovar.coding_interval_support import AlleleReadTiling, coding_interval_support
from isovar.dna import reverse_complement_dna
from isovar.rna_evidence import EvidenceSets


@pytest.mark.parametrize("strand", ["+", "-"])
def test_allele_anchored_partial_single_and_merged_support(strand):
    target = "ATGCAAGGCTAA"
    observations = [
        ("full", "ATG", "C", "AAGGCTAA", 1),
        ("partial", "G", "C", "AAG", 1),
        ("pair", "ATG", "C", "AAGGCTAA", 2),
        ("wrong", "ATG", "C", "AAGACTAA", 1),
    ]
    reads = []
    for name, prefix, allele, suffix, count in observations:
        if strand == "-":
            prefix, allele, suffix = map(reverse_complement_dna, (suffix, allele, prefix))
        reads.append(AlleleRead(prefix, allele, suffix, name, source_read_count=count))
    evidence = EvidenceSets(["sample", "library"])
    result = coding_interval_support(reads + reads[:1], target, 3, (0, 12), strand, evidence.support)
    assert result["allele_observation_span_support"]["fragments"] == 2
    assert result["single_read_span_support"]["fragments"] == 1
    assert result["paired_template_span_support"]["fragments"] == 1
    assert result["fragment_depth"] == [2, 2, 3, 3, 3, 3, 3, 2, 2, 2, 2, 2]
    assert result["minimum_fragment_depth"] == 2
    assert result["fully_tiled"]


def test_empty_interval_does_not_claim_complete_sequence():
    result = coding_interval_support([], "", 0, (0, 0), "+", EvidenceSets(["s", "r"]).support)
    assert result["fragment_depth"] == []
    assert not result["fully_tiled"]


def test_single_and_merged_subsets_can_share_a_fragment():
    one = AlleleRead("ATG", "C", "AAG", "template", source_read_count=1)
    pair = AlleleRead("ATG", "C", "AAG", "template", source_read_count=2)
    result = coding_interval_support([one, pair], "ATGCAAG", 3, (0, 7), "+", EvidenceSets(["s", "r"]).support)
    assert result["allele_observation_span_support"]["fragments"] == 1
    assert result["single_read_span_support"]["fragments"] == result["paired_template_span_support"]["fragments"] == 1


def test_mismatch_excludes_a_read_only_from_intervals_that_contain_it():
    target = "ATGCAAGGCTAA"
    clean = AlleleRead("ATG", "C", "AA", "template")
    mismatched = AlleleRead("ATG", "C", "AAGGCTAC", "template")  # disagrees at the last base
    other = AlleleRead("ATG", "C", "AAGGCTAC", "other")
    tiling = AlleleReadTiling([clean, mismatched, other], target, 3, "+")
    support = EvidenceSets(["s", "r"]).support
    inside = tiling.interval_support((0, 9), support)
    assert inside["fragment_depth"] == [2] * 9
    assert inside["allele_observation_span_support"]["reads"] == 2
    whole = tiling.interval_support((0, 12), support)
    assert whole["fragment_depth"] == [1] * 6 + [0] * 6
    assert whole["allele_observation_span_support"]["reads"] == 0
    assert not whole["fully_tiled"]


def _depths_by_definition(observations, target, variant_start, start, end):
    """Fragments with a read that covers a base and agrees with the target inside the interval."""
    covering = [set() for _ in range(start, end)]
    for name, prefix, allele, suffix, _ in observations:
        text, first = prefix + allele + suffix, variant_start - len(prefix)
        inside = [offset for offset in range(start, end) if first <= offset < first + len(text)]
        if all(text[offset - first] == target[offset] for offset in inside):
            for offset in inside:
                covering[offset - start].add(name)
    return [len(names) for names in covering]


@pytest.mark.parametrize("strand", ["+", "-"])
def test_one_tiling_answers_every_interval_like_a_fresh_measurement(strand):
    target = "ATGCAAGGCTAA"
    observations = [("a", "ATG", "C", "AAGGCTAA", 1), ("a", "G", "C", "AAGTC", 2), ("b", "TG", "C", "AAGGC", 1),
                    ("c", "", "C", "AAGGCTAAGG", 1), ("b", "TTATG", "C", "A", 2), ("legacy", "ATG", "C", "TAG", 1)]
    reads = []
    for index, (name, prefix, allele, suffix, count) in enumerate(observations):
        if strand == "-":
            prefix, allele, suffix = map(reverse_complement_dna, (suffix, allele, prefix))
        alignments = () if name == "legacy" else [(("group", name, 0), (0, index, "8M", False))]
        reads.append(AlleleRead(prefix, allele, suffix, name, source_read_count=count, source_alignments=alignments))
    tiling = AlleleReadTiling(reads, target, 3, strand)
    for start in range(len(target) + 1):
        for end in range(start, len(target) + 1):
            shared, fresh = EvidenceSets(["s", "r"]), EvidenceSets(["s", "r"])
            result = tiling.interval_support((start, end), shared.support)
            assert result == coding_interval_support(reads, target, 3, (start, end), strand, fresh.support)
            assert result["fragment_depth"] == _depths_by_definition(observations, target, 3, start, end)


@pytest.mark.parametrize("arguments", [("ATGC", 5, "+"), ("ATGC", 1, "?")])
def test_invalid_tiling_is_rejected(arguments):
    with pytest.raises(ValueError):
        AlleleReadTiling([], *arguments)
    with pytest.raises(ValueError):
        AlleleReadTiling([], "ATGC", 1, "+").interval_support((3, 2), EvidenceSets(["s", "r"]).support)
