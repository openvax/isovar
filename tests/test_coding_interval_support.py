"""Partial contributors and merged mates are not single-read coding spans."""

import pytest

from isovar.allele_read import AlleleRead
from isovar.coding_interval_support import coding_interval_support
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
