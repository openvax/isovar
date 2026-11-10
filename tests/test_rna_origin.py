"""Competitive origin is molecule-linked context, never a contig/code shortcut."""

from dataclasses import replace
import json
import random

import pysam
import pytest

from isovar import RNAOriginHypothesis, RNAOriginPolicy, assess_rna_origin, prepare_rna_library
from isovar.dna import reverse_complement_dna
from isovar.read_collector import ReadCollector
from isovar.rna_path_geometry import reverse_path
from tests.test_sv_rna import aligned, record, write_bam


def change(base):
    return next(b for b in "ACGT" if b != base)


@pytest.fixture
def panel():
    rng = random.Random(704)
    mt = "".join(rng.choice("ACGT") for _ in range(200))
    nuc = list(mt)
    for i in (40, 60, 140, 160):
        nuc[i] = change(nuc[i])
    return [RNAOriginHypothesis("mt", "mitochondrial", mt,
                               [("1", p, "+") for p in range(1000, 1200)], (100, 101)),
            RNAOriginHypothesis("numt", "nuclear", "".join(nuc),
                               [("2", p, "+") for p in range(2000, 2200)], (100, 101))]


def observed(hypothesis, name="private-original", *, start=20, end=180, flag=0, group="a", **kwargs):
    # Deliberately map both sequence origins to the nominated mitochondrial locus.
    return record(name, "1", 1000 + start, f"{end - start}M", hypothesis.sequence[start:end],
                  flag=flag, group=group, **kwargs)


def assess(tmp_path, panel, reads, *, source="L", policy=None):
    path = write_bam(tmp_path / "origin.bam", reads)
    with pysam.AlignmentFile(path) as bam:
        return assess_rna_origin(bam, panel, sample_id="S", source=source, policy=policy)


@pytest.mark.parametrize("origin,index", [("mitochondrial", 0), ("nuclear", 1)])
@pytest.mark.parametrize("flag", [0, 16])
@pytest.mark.parametrize("mapq", [0, 255])
def test_sequence_origin_independent_of_mapping_strand_and_mapq(tmp_path, panel, origin, index, flag, mapq):
    result = assess(tmp_path, panel, [observed(panel[index], flag=flag, mapq=mapq)])
    assert result["assessment"] == origin + "_favored"
    assert result["rna_support"]["reads"] == result["rna_support"]["fragments"] == 1
    assert result["support_by_status"][origin + "_favored"]["fragments"] == 1
    saved = json.loads(json.dumps(result))
    assert saved["panel_id"] == result["panel_id"]
    assert "private-original" not in json.dumps(saved)
    assert "translation_table" not in saved
    for evidence in result["evidence_sets"].values():
        assert len(set(evidence["read_ids"])) == evidence["reads"]
        assert len(set(evidence["fragment_ids"])) == evidence["fragments"]


def test_reverse_oriented_nuclear_hypothesis(tmp_path, panel):
    sequence, positions = reverse_path(panel[1].sequence, panel[1].positions)
    reverse = replace(panel[1], sequence=sequence, positions=positions, target_interval=(99, 100))
    result = assess(tmp_path, [panel[0], reverse], [observed(panel[1])])
    assert result["assessment"] == "nuclear_favored"


def test_focal_allele_is_not_an_origin_marker(tmp_path, panel):
    # Neither the original base nor its reverse complement.
    alternative = "C" if panel[0].sequence[100] in "AT" else "A"
    changed = [replace(h, sequence=h.sequence[:100] + alternative + h.sequence[101:]) for h in panel]
    result = assess(tmp_path, changed, [observed(changed[1])])
    comparison = result["fragments"][0]["segments"][0]["observations"][0]
    assert result["assessment"] == "nuclear_favored"
    assert comparison["origin_edit_distances"]["nuclear"] == 0
    assert comparison["edit_margin"] >= 3
    with pytest.raises(ValueError, match="same target"):
        assess(tmp_path, [panel[0], changed[1]], [])


def test_identical_competitors_stay_ambiguous(tmp_path, panel):
    panel[1] = replace(panel[1], sequence=panel[0].sequence)
    assert assess(tmp_path, panel, [observed(panel[0])])["assessment"] == "ambiguous"


@pytest.mark.parametrize("limitation", ["short", "missing_qualities", "low_target", "poor_fit", "short_panel"])
def test_unassessed_is_not_negative_evidence(tmp_path, panel, limitation):
    read = observed(panel[0], qualities=limitation != "missing_qualities")
    if limitation == "short":
        read = observed(panel[0], start=90, end=111)
    elif limitation == "low_target":
        qualities = list(read.query_qualities)
        qualities[80] = 1
        read.query_qualities = qualities
    elif limitation == "poor_fit":
        sequence = list(read.query_sequence)
        for i in range(0, 60, 3):
            sequence[i] = change(sequence[i])
        read.query_sequence = "".join(sequence)
        read.query_qualities = [30] * len(sequence)
    elif limitation == "short_panel":
        panel[1] = replace(panel[1], sequence=panel[1].sequence[50:151],
                           positions=panel[1].positions[50:151], target_interval=(50, 51))
    result = assess(tmp_path, panel, [read])
    assert result["assessment"] == "unassessed"
    assert result["support_by_status"]["unassessed"]["fragments"] == 1
    assert result["rna_support"]["fragments"] == 1


def test_low_quality_context_cannot_supply_discriminating_markers(tmp_path, panel):
    read = observed(panel[1])
    qualities = list(read.query_qualities)
    for i in (40, 60, 140, 160):
        qualities[i - 20] = 1
    read.query_qualities = qualities
    result = assess(tmp_path, panel, [read])
    assert result["assessment"] == "ambiguous"
    scores = result["fragments"][0]["segments"][0]["observations"][0]
    assert scores["masked_context_bases"] == 4
    assert scores["required_edit_margin"] == 7


def test_original_records_placements_mates_and_groups_are_not_votes(tmp_path, panel):
    read = observed(panel[1], name="repeated")
    reads = [read, read, observed(panel[1], name="repeated", flag=256)]
    reads += [observed(panel[1], name="pair", flag=65), observed(panel[1], name="pair", flag=129)]
    reads += [observed(panel[1], name="repeated", group="b")]
    result = assess(tmp_path, panel, reads)
    assert result["assessment"] == "nuclear_favored"
    assert result["rna_support"]["reads"] == 4
    assert result["rna_support"]["fragments"] == 3
    assert len(result["fragments"]) == 3
    assert sum(row["reads"] for row in result["support_by_status"].values()) == 4


@pytest.mark.parametrize("same_segment", [False, True])
def test_conflicting_mates_or_original_placements_are_retained(tmp_path, panel, same_segment):
    reads = [observed(panel[0], flag=0 if same_segment else 65),
             observed(panel[1], flag=256 if same_segment else 129)]
    result = assess(tmp_path, panel, reads)
    assert result["assessment"] == "conflicting"
    assert result["rna_support"]["fragments"] == 1
    assert result["rna_support"]["reads"] == (1 if same_segment else 2)


def test_disjoint_mate_cannot_create_unobserved_allele_phase(tmp_path, panel):
    reads = [observed(panel[0], flag=65, start=80, end=121),
             observed(panel[1], flag=129, start=20, end=80)]
    result = assess(tmp_path, panel, reads)
    assert result["assessment"] == "ambiguous"
    assert result["rna_support"]["reads"] == 1


def test_mixed_origins_are_not_pooled_into_one_call(tmp_path, panel):
    reads = [observed(panel[0], name="mt"), observed(panel[1], name="nuc")]
    result = assess(tmp_path, panel, reads)
    assert result["assessment"] == "mixed_origins"
    assert {r["status"] for r in result["fragments"]} == {"nuclear_favored", "mitochondrial_favored"}
    other_library = assess(tmp_path, panel, reads, source="L2")
    assert {f["fragment_id"] for f in result["fragments"]}.isdisjoint(
        f["fragment_id"] for f in other_library["fragments"])


@pytest.mark.parametrize("uninformative,status", [
    # Shared sequence only: no distinguishing marker is reached.
    (dict(start=70, end=130), "ambiguous"),
    (dict(start=70, end=130, qualities=None), "unassessed"),
])
def test_uninformative_fragments_are_counted_but_cannot_veto_the_sample(tmp_path, panel, uninformative, status):
    reads = [observed(panel[1], name=f"nuc{i}") for i in range(3)]
    reads.append(observed(panel[1], name="shared", **uninformative))
    result = assess(tmp_path, panel, reads)
    assert result["assessment"] == "nuclear_favored"
    counts = {s: row["fragments"] for s, row in result["support_by_status"].items()}
    assert counts["nuclear_favored"] == 3 and counts[status] == 1
    assert assess(tmp_path, panel, [reads[-1]])["assessment"] == status
    conflicted = reads + [observed(panel[0], name="pair", flag=65), observed(panel[1], name="pair", flag=129)]
    assert assess(tmp_path, panel, conflicted)["assessment"] == "conflicting"


def test_wrong_allele_does_not_enter_the_origin_denominator(tmp_path, panel):
    wrong = replace(panel[0], sequence=panel[0].sequence[:100] + change(panel[0].sequence[100]) + panel[0].sequence[101:])
    result = assess(tmp_path, panel, [observed(wrong)])
    assert result["assessment"] == "no_target_support"
    assert result["rna_support"]["reads"] == 0


def test_competing_allele_placement_cannot_be_discarded(tmp_path, panel):
    wrong = replace(panel[0], sequence=panel[0].sequence[:100] + change(panel[0].sequence[100]) + panel[0].sequence[101:])
    result = assess(tmp_path, panel, [observed(panel[1]), observed(wrong, flag=256)])
    assert result["assessment"] == "ambiguous"


@pytest.mark.parametrize("gap", [3, 1000])
def test_anchored_deletion_or_splice_target_uses_observed_path(tmp_path, panel, gap):
    panel = [replace(h, positions=[(p[0], p[1] + (gap if i >= 100 else 0), p[2])
                                   for i, p in enumerate(h.positions)], target_interval=(95, 105)) for h in panel]
    reads = aligned("junction", panel[1].sequence[20:180], panel[1].positions[20:180])
    result = assess(tmp_path, panel, reads)
    assert result["assessment"] == "nuclear_favored"
    # CIGAR N cannot manufacture observed target bases inside a skipped region.
    skipped = record("skipped", "1", 1020, "60M1000N100M", panel[1].sequence[20:180])
    assert assess(tmp_path, panel, [skipped])["assessment"] == "no_target_support"


def test_nuclear_exon_context_and_reused_library_index(tmp_path, panel):
    rng = random.Random(100)
    nuclear = panel[0].sequence[:120] + "".join(rng.choice("ACGT") for _ in range(80))
    panel[1] = replace(panel[1], sequence=nuclear, positions=[
        ("2", p + (1000 if i >= 120 else 0), "+") for i, p in enumerate(range(2000, 2200))])
    path = write_bam(tmp_path / "nuclear.bam", aligned("spliced", nuclear[20:180], panel[1].positions[20:180]))
    prepared = prepare_rna_library(path, tmp_path / "index", read_collector=ReadCollector.for_rna_paths())
    with pysam.AlignmentFile(path) as bam, prepared.open() as index:
        result = assess_rna_origin(bam, panel, sample_id="S", source="L", library_index=index)
    assert result["assessment"] == "nuclear_favored"
    assert result["provenance"]["input_identity"] == prepared.identity


def test_anchored_insertion_target(tmp_path, panel):
    panel = [replace(h, sequence=h.sequence[:100] + "ACG" + h.sequence[100:],
                     positions=h.positions[:100] + (None,) * 3 + h.positions[100:],
                     target_interval=(99, 104)) for h in panel]
    read = record("inserted", "1", 1020, "80M3I80M", panel[1].sequence[20:183])
    result = assess(tmp_path, panel, [read])
    assert result["assessment"] == "nuclear_favored"
    assert result["rna_support"]["fragments"] == 1


@pytest.mark.parametrize("linkage", ["actual", "SA_only", "unlinked"])
def test_nuclear_flank_requires_actual_linked_supplementary_records(tmp_path, panel, linkage):
    rng = random.Random(101)
    nuclear = panel[0].sequence[:120] + "".join(rng.choice("ACGT") for _ in range(80))
    # The observed target lies in the shared mitochondrial sequence; only
    # the linked distal nuclear flank distinguishes its origin.
    panel[1] = replace(panel[1], sequence=nuclear,
                       positions=panel[0].positions[:120] + panel[1].positions[120:])
    reads = aligned("split", nuclear[20:180], panel[1].positions[20:180])
    if linkage == "SA_only":
        reads = [r for r in reads if not r.is_supplementary]
    elif linkage == "unlinked":
        for read in reads:
            read.set_tag("SA", None)
    result = assess(tmp_path, panel, reads)
    assert result["assessment"] == ("nuclear_favored" if linkage == "actual" else "ambiguous")
    assert result["rna_support"]["reads"] == result["rna_support"]["fragments"] == 1


def test_competitor_order_and_duplicates_cannot_outvote_a_tie(tmp_path, panel):
    reads = [observed(panel[0])]
    panel += [replace(panel[0], identifier="other-numt", origin="nuclear")]
    result = assess(tmp_path, panel, reads)
    assert result["assessment"] == "ambiguous"
    assert assess(tmp_path, panel[::-1], reads) == result


def test_complementary_targets_at_same_placement_remain_alternatives(tmp_path, panel):
    target = reverse_complement_dna(panel[0].target_sequence)
    panel[1] = replace(panel[1], positions=panel[0].positions,
                       sequence=panel[1].sequence[:100] + target + panel[1].sequence[101:])
    reads = [observed(panel[0])]
    assert assess(tmp_path, panel, reads)["assessment"] == "ambiguous"
    renamed = [replace(panel[0], identifier="z"), replace(panel[1], identifier="a")]
    assert assess(tmp_path, renamed, reads)["assessment"] == "ambiguous"


def test_cli_roundtrip_and_existing_output_protection(tmp_path, panel):
    from isovar.cli.commands import run

    bam = write_bam(tmp_path / "input.bam", [observed(panel[1])])
    request, output = tmp_path / "panel.json", tmp_path / "result.json"
    request.write_text(json.dumps(dict(schema="isovar.rna_origin_panel.v1", hypotheses=[h.to_dict() for h in panel])))
    args = ["rna-origin", "--bam", str(bam), "--hypotheses", str(request), "--sample-id", "S",
            "--source-id", "L", "--output", str(output)]
    run(args)
    result = json.loads(output.read_text())
    assert result["schema"] == "isovar.rna_origin.v1"
    assert result["assessment"] == "nuclear_favored"
    with pytest.raises(SystemExit):
        run(args)


@pytest.mark.parametrize("kwargs", [{"min_edit_margin": 0}, {"min_context_bases": True},
    {"min_base_quality": -1}, {"max_edit_fraction": float("nan")}, {"max_edit_fraction": 1}])
def test_invalid_policy(kwargs):
    with pytest.raises(ValueError):
        RNAOriginPolicy(**kwargs)


def test_invalid_or_incomplete_panel(tmp_path, panel):
    with pytest.raises(ValueError, match="RNAOriginHypothesis"):
        assess(tmp_path, [None], [])
    with pytest.raises(ValueError, match="RNAOriginPolicy"):
        assess(tmp_path, panel, [], policy={})
    with pytest.raises(ValueError, match="nuclear and mitochondrial"):
        assess(tmp_path, panel[:1], [])
    with pytest.raises(ValueError, match="unique"):
        assess(tmp_path, [panel[0], panel[0], panel[1]], [])
    with pytest.raises(ValueError, match="placed anchor"):
        replace(panel[0], positions=[None] * len(panel[0].sequence))
    with pytest.raises(ValueError, match="outside"):
        assess(tmp_path, [panel[0], replace(panel[1], positions=[("missing", i, "+") for i in range(200)])], [])
