"""The catalogue audit must preserve its denominator and evidence semantics."""

from copy import deepcopy
import gzip
import json
from pathlib import Path

import osteosarc
import pytest
import pysam

from examples.sid_sv_audit.inventory import (
    classify_source, fusion_breakends, geometry, nominations, oriented_breakpoints,
    load_inventory,
)
from examples.sid_sv_audit.report import coverage
from examples.sid_sv_audit.reconstruct import reconstruct_one, summarize
from isovar import sid_data


def test_retained_flanks_determine_both_views_including_inversions():
    ends = fusion_breakends(dict(chrD="chr2", posD="101", retainD="right",
                                chrA="chr7", posA="201", retainA="right"))
    assert ends == [dict(contig="chr2", position=100, retained_side="right"),
                    dict(contig="chr7", position=200, retained_side="right")]
    assert oriented_breakpoints(ends) == (dict(contig="chr2", position=100, strand="-"),
                                         dict(contig="chr7", position=200, strand="+"))
    assert oriented_breakpoints(ends, reverse=True) == (dict(contig="chr7", position=200, strand="-"),
                                                       dict(contig="chr2", position=100, strand="+"))
    assert geometry(ends[::-1]) == ends
    with pytest.raises(ValueError, match="Invalid"):
        geometry([dict(ends[0], retained_side=None), ends[1]])


def test_entire_shared_catalogue_retains_unresolved_targets_and_aliases():
    catalogue = osteosarc.load_sv_candidates()
    targets, groups, survey = nominations(catalogue, [], [])
    assert len(targets) == 637
    assert survey == []
    assert sum(t["geometry_id"] is None for t in targets.values()) == 120
    assert {name for group in groups.values() for name in group["nominations"]} == {
        name for name, target in targets.items() if target["geometry_id"] is not None}
    assert targets["shared:SV0001"]["original"]["original_calls"][0]["vcf"]["alt"].startswith("GTGTGCTGATATG[")
    # PARD3B's acceptor begins before the last retained base, not after it.
    pard = groups[targets["shared:SV0001"]["geometry_id"]]
    assert pard["breakends"][1] == dict(contig="chr9", position=22007647, retained_side="right")


def test_frozen_full_catalogue_and_source_denominator():
    directory = Path(__file__).parents[1] / "examples" / "sid_sv_audit" / "data"
    manifest = load_inventory(directory)
    assert len(manifest["nominations"]) == 2258
    assert len(manifest["fusion_survey"]) == 1271
    assert len(manifest["geometries"]) == 1493
    assert len(manifest["sources"]) == 913
    assert sum(s["selection"]["cohort"] == "tumor_candidate" for s in manifest["sources"].values()) == 77
    assert {n for g in manifest["geometries"].values() for n in g["nominations"]} == {
        n for n, t in manifest["nominations"].items() if t["geometry_id"]}
    sources = {n: s for n, s in manifest["sources"].items() if s["selection"]["cohort"] == "tumor_candidate"}
    ledger = coverage(manifest, sources, {}, {})
    assert ledger["expected_pairs"] == 173866
    assert len(ledger["rows"]) == ledger["expected_pairs"]
    assert not ledger["all_pairs_accounted"]


def test_survey_calls_are_not_filtered_by_support_or_collapsed_into_alleles():
    row = dict(name="A::B", caller="arriba", sample="T0", chrD="chr2", posD="101", retainD="left",
               chrA="chr7", posA="201", retainA="right", insert="A")
    second = dict(row, caller="linx_sv", insert="G")
    summary = [dict(locusA="chr2:101", locusB="chr7:201", label="A::B", verdict="NONE")]
    targets, groups, _ = nominations(dict(targets={}), [row, second], summary)
    assert len(targets) == 2 and len(groups) == 1
    assert {t["original"]["insert"] for t in targets.values()} == {"A", "G"}
    with pytest.raises(ValueError, match="inventory disagree"):
        nominations(dict(targets={}), [row], [])


def test_source_scope_keeps_blood_separate_and_does_not_turn_a_path_into_metadata():
    row = dict(key="rna-seq/reprocessed/unknown/sample.bam", claims=[])
    selection = classify_source(row)
    assert selection["cohort"] == "tumor_candidate"
    assert selection["library_claims"] == selection["timepoint_claims"] == []
    assert selection["reason"] == "RNA_path_hint"
    assert classify_source(dict(key="rna-seq/tempus/name/WES/sample.bam", claims=[]))["cohort"] == "excluded"
    assert classify_source(dict(key="kamil/tumor/output/T2/outs/vdj_b/consensus.bam", claims=[
        dict(assay="scrna-seq", tissue="tumor", basis="published")]))["reason"] == "targeted_immune_repertoire_product"
    blood = dict(key="kamil/blood/output/Jan2025_GEX/outs/possorted_genome_bam.bam", claims=[
        dict(assay="scrna-seq", tissue="blood", basis="published")])
    assert classify_source(blood)["cohort"] == "blood_control"


def test_coverage_distinguishes_missing_work_from_unresolved_and_unavailable():
    manifest = dict(nominations={"a": dict(geometry_id="g"), "alias": dict(geometry_id="g"),
                                 "single": dict(geometry_id=None)}, geometries={"g": {}})
    sources = {"one": {}, "two": {}}
    headers = {"one": dict(status="ready"), "two": dict(status="no_listed_index")}
    outcomes = {("one", "g", "forward"): dict(status="no_candidate_paths")}
    result = coverage(manifest, sources, outcomes, headers)
    assert result["expected_pairs"] == 6
    assert not result["all_pairs_accounted"]
    assert result["counts"] == dict(pending_reconstruction=2, unresolved_geometry=2, not_assessable=2)
    outcomes["one", "g", "reverse"] = dict(status="event_linked_candidates")
    result = coverage(manifest, sources, outcomes, headers)
    assert result["all_pairs_accounted"]
    assert result["observed_geometry_views"] == 2  # Never count the alias again.
    outcomes["extraneous", "g", "forward"] = dict(status="no_candidate_paths")
    with pytest.raises(ValueError, match="outside"):
        coverage(manifest, sources, outcomes, headers)


def test_new_audit_reconstructs_original_pacbio_orf_without_qualities(tmp_path):
    """Independent pinned protein/count expectation, through the new adapter."""
    root = Path(__file__).parent / "data" / "fusions"
    with gzip.open(root / "corpus" / "TPST1--CRCP-T1.input.json.gz", "rt") as handle:
        original = json.load(handle)
    fusion = original["fusion"]
    ends = [dict(contig=fusion["donor"]["contig"], position=fusion["donor"]["position"],
                 retained_side="left" if fusion["donor"]["strand"] == "+" else "right"),
            dict(contig=fusion["acceptor"]["contig"], position=fusion["acceptor"]["position"],
                 retained_side="right" if fusion["acceptor"]["strand"] == "+" else "left")]
    group = dict(breakends=ends, nominations=["TPST1--CRCP"])
    reference = dict(models={r["transcript_id"]: r for r in original["references"]},
                     assignments={"test": dict(transcripts=[r["transcript_id"] for r in original["references"]])})
    source = dict(id="pacbio-fixture", url="pinned-original-PacBio-T1", claims=[])
    bam_path = tmp_path / "original.bam"
    sid_data.export("fusions/long-read/TPST1--CRCP.PacBio-T1.sam.gz", bam_path)
    with pysam.AlignmentFile(str(bam_path)) as bam:
        result = reconstruct_one(bam, "test", group, source, reference, "forward", dict(assemble=False))
    summary = summarize(result)
    candidate, = [c for c in summary["candidates"] if c["amino_acids"] == "MPSRRRGGSSLIPGRMPILRFSVTTRYFSY"]
    assert candidate["ends_with_stop_codon"]
    assert candidate["rna_support"]["fragments"] == 15
    assert candidate["rna_support"]["missing_quality_reads"] == 15
    assert not candidate["initiation_observed"] and not candidate["translation_observed"]
    # Regional copies of an identical ORF must not inflate event attribution.
    before = deepcopy(summary)
    duplicate = deepcopy(result["paths"][0])
    for junction in duplicate["junctions"]:
        junction["relation"] = "regional_novel_junction"
    result["paths"].append(duplicate)
    assert summarize(result)["candidates"] == before["candidates"]
