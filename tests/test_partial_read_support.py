"""Conserved fragment budgets distinguish informative and ambiguous partial reads."""

from copy import copy, deepcopy
import json
from math import isclose
from types import SimpleNamespace

import pytest

from isovar import AlleleRead, BaseQualityPolicy
from isovar.base_quality import AlleleQualityFootprint
from isovar.protein_hypotheses import _EvidenceSets
from isovar.rna_candidate_support import score_fragments

TX = SimpleNamespace(id="tx", exon_intervals=[(1, 20000)])


def candidate(key, suffix="TTTAAA", protein=None, prefix="AAACCC"):
    return dict(candidate_id=key, sequence=prefix + "G" + suffix,
                variant_interval=[len(prefix), len(prefix) + 1],
                transcript_ids=["tx"], hypothesis_ids=[protein or key])


def read(name, prefix="CCC", allele="G", suffix="TTT", mate=0, group="rg", placement=100, quality=None):
    identity = ((group, name, mate), (0, placement, "%dM" % (len(prefix) + len(allele) + len(suffix)), False))
    footprint = AlleleQualityFootprint(identity, allele, (len(prefix), len(prefix) + len(allele)),
                                       (len(prefix), len(prefix) + len(allele)), (len(prefix),),
                                       (quality,), "allele_bases")
    return AlleleRead(prefix, allele, suffix, name, source_alignments=(identity,),
                     source_allele_qualities=(footprint,))


def score(reads, candidates=None, **kwargs):
    candidates = candidates if candidates is not None else [candidate("a"), candidate("b", suffix="TTTCCC")]
    policy = kwargs.pop("quality_policy", None)
    evidence = _EvidenceSets(["s", "input"], policy)
    result = score_fragments(reads, candidates, {c["candidate_id"]: [TX] for c in candidates}, "A", evidence,
                             min_overlap=3, **kwargs)
    assert sum(s["fractional_fragments"] for s in result["protein_scores"]) == pytest.approx(result["scored_fragments"])
    assert sum(s["fractional_fragments"] for s in result["candidate_scores"]) == pytest.approx(result["scored_fragments"])
    return result


def test_short_informative_read_keeps_full_weight_and_shared_reads_split_it():
    reads = [read("shared%d" % i) for i in range(9)] + [read("informative", suffix="TTTAAA")]
    result = score(reads)
    assert result["input_fragments"] == result["scored_fragments"] == 10
    a, b = result["candidate_scores"]
    assert (a["fractional_fragments"], a["unique_fragments"], a["ambiguous_fragments"], a["rank"]) == (5.5, 1, 9, 1)
    assert (b["fractional_fragments"], b["unique_fragments"], b["ambiguous_fragments"], b["rank"]) == (4.5, 0, 9, 2)
    informative = next(f for f in result["fragments"] if f["candidate_ids"] == ["a"])
    assert informative["candidate_weight"] == 1
    assert not informative["observations"][0]["comparisons"][0][0]["spans_candidate"]
    assert json.loads(json.dumps(result)) == result
    assert result["policy"]["calibrated_probability"] is False


def test_unobserved_extensions_are_neutral_and_near_matches_remain_alternatives():
    reads = [read("short")]
    original = [candidate("a"), candidate("b", suffix="TCTCCC")]
    # Both candidates remain compatible at tolerance 1; minimum-cost selection
    # would turn a noisy single base into unjustified unique support.
    initial = score(reads, original)
    assert [s["fractional_fragments"] for s in initial["candidate_scores"]] == [.5, .5]
    extended = deepcopy(original)
    for c in extended:
        c["sequence"] = "T" * 50 + c["sequence"] + "C" * 50
        c["variant_interval"] = [i + 50 for i in c["variant_interval"]]
    assert score(reads, extended)["candidate_scores"] == initial["candidate_scores"]
    exact = score(reads, original, max_edits=0)
    assert [s["fractional_fragments"] for s in exact["candidate_scores"]] == [1, 0]


def test_fragments_not_mates_records_or_read_group_name_collisions():
    first = read("pair", suffix="TTTAAA", mate=64)
    second = read("pair", mate=128)
    separate = score([first, second, first, second])
    assert separate["input_fragments"] == separate["scored_fragments"] == 1
    assert separate["fragments"][0]["rna_support"]["reads"] == 2
    merged = AlleleRead(first.prefix, "G", first.suffix, "pair", source_read_count=2,
                        source_alignments=first.source_alignments + second.source_alignments)
    combined = score([merged, first, second])
    assert combined["candidate_scores"] == separate["candidate_scores"]
    two_groups = score([first, second, read("pair", group="different")])
    assert two_groups["scored_fragments"] == 2
    assert len({f["fragment_id"] for f in two_groups["fragments"]}) == 2


def test_mate_conflicts_and_disagreeing_placements_are_explicitly_unscored():
    conflict = score([read("pair", suffix="TTTAAA", mate=64), read("pair", suffix="TTTCCC", mate=128)])
    assert conflict["unscored_fragments"] == 1
    assert conflict["fragments"][0]["status"] == "conflicting_observations"
    alternative = score([read("one", suffix="TTTAAA"), read("one", suffix="TTTCCC", placement=200)])
    assert alternative["input_fragments"] == 1
    assert alternative["fragments"][0]["status"] == "alternative_placements_disagree"
    agreeing = score([read("one"), read("one", placement=200)])
    assert agreeing["scored_fragments"] == 1
    assert [s["fractional_fragments"] for s in agreeing["candidate_scores"]] == [.5, .5]


def test_reference_competitors_short_reads_and_empty_catalogs_are_accounted_for():
    reads = [read("ref", allele="A"), read("other", allele="T"), read("short", prefix="", suffix=""),
             read("unknown", prefix="NNN"), read("alt")]
    result = score(reads)
    assert (result["input_fragments"], result["scored_fragments"], result["unscored_fragments"]) == (5, 1, 4)
    assert all(f["candidate_weight"] == f["protein_weight"] == 0 for f in result["fragments"] if f["status"] != "scored")
    empty = score(reads, [])
    assert empty["input_fragments"] == empty["unscored_fragments"] == 5
    assert {f["status"] for f in empty["fragments"]} == {"no_candidate"}
    assert score([])["input_fragments"] == 0


def test_synonymous_rna_candidates_do_not_dilute_or_inflate_protein_weights():
    candidates = [candidate("a", protein="P"), candidate("b", protein="Q", suffix="TTTCCC")]
    original = score([read("short")], candidates)
    duplicate_protein = candidate("synonymous", protein="P", prefix="GGGCCC")
    expanded = score([read("short")], candidates + [duplicate_protein])
    assert expanded["protein_scores"] == original["protein_scores"]
    assert all(isclose(s["fractional_fragments"], 1 / 3) for s in expanded["candidate_scores"])
    only_one_protein = score([read("short")], [candidates[0], duplicate_protein])
    assert only_one_protein["protein_scores"][0]["unique_fragments"] == 1
    assert only_one_protein["protein_scores"][0]["fractional_fragments"] == 1


def test_quality_strata_preserve_score_and_unknowns_without_platform_defaults():
    reads = [read("high", quality=40), read("low", quality=0), read("missing"),
             read("pair", mate=64, quality=40), read("pair", mate=128)]
    raw = score(reads)
    qualified = score(reads, quality_policy=BaseQualityPolicy(20))
    for score_row, raw_score in zip(qualified["protein_scores"], raw["protein_scores"]):
        assert {k: v for k, v in score_row.items() if k != "base_quality_fractional_fragments"} == raw_score
        assert score_row["base_quality_fractional_fragments"] == {"passed": .5, "failed": .5, "unassessed": 1}
    # Stripping all quality provenance leaves every assignment and rank intact.
    missing = [copy(r) for r in reads]
    for r in missing:
        r.source_allele_qualities = ()
    assert score(missing)["protein_scores"] == raw["protein_scores"]


def test_deterministic_scores_dense_ties_and_incomplete_catalog_metadata():
    reads = [read(str(i)) for i in range(7)]
    candidates = [candidate("a"), candidate("b"), candidate("c")]
    result = score(reads, candidates, catalog_complete=False)
    assert result["catalog_complete"] is False
    assert {s["rank"] for s in result["protein_scores"]} == {1}
    assert score(list(reversed(reads)), list(reversed(candidates)))["fragments"] == result["fragments"]
    assert score(reads + reads, candidates)["protein_scores"] == result["protein_scores"]


def test_legacy_reads_keep_null_identity_and_quality():
    old = AlleleRead("CCC", "G", "TTT", "legacy")
    result = score([old], quality_policy=BaseQualityPolicy(20))
    fragment, = result["fragments"]
    assert fragment["fragment_id"] is fragment["rna_support"]["evidence_set_id"] is None
    assert fragment["base_quality_status"] == "unassessed"
    assert result["scored_fragments"] == 1


@pytest.mark.parametrize("kwargs", [dict(max_edits=-1), dict(max_edits=True)])
def test_invalid_comparison_settings(kwargs):
    with pytest.raises(ValueError):
        score([], **kwargs)
