"""The catalogue audit must preserve its denominator and evidence semantics."""

from copy import deepcopy
import gzip
import json
from pathlib import Path
import shutil
from types import SimpleNamespace

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

from examples.sid_sv_audit import acquisition, inventory, reconstruct, report
from examples.sid_sv_audit.inventory import digest, identity, read_json, write_json


def fake_subset(bam, receipt):
    path = bam.parent / "receipt.json"
    write_json(path, receipt)
    return osteosarc.ReadSubset(bam, Path(str(bam) + ".bai"), receipt)


@pytest.mark.parametrize("name", ["result.json", "result.json.gz"])
def test_checkpoint_write_interruption_never_publishes_partial_content(tmp_path, monkeypatch, name):
    path = tmp_path / name
    original = inventory.tempfile.NamedTemporaryFile

    class InterruptedWrite:
        def __init__(self, **kwargs):
            self.handle = original(**kwargs)
            self.name = self.handle.name

        def __enter__(self):
            return self

        def __exit__(self, *args):
            self.handle.close()

        def write(self, data):
            self.handle.write(data[:len(data) // 2])
            self.handle.flush()
            raise KeyboardInterrupt()

    monkeypatch.setattr(inventory.tempfile, "NamedTemporaryFile", InterruptedWrite)
    with pytest.raises(KeyboardInterrupt):
        write_json(path, {"complete": True})
    assert not path.exists()
    assert not list(tmp_path.glob(".checkpoint-*"))
    monkeypatch.setattr(inventory.tempfile, "NamedTemporaryFile", original)
    write_json(path, {"complete": True})
    write_json(path, {"complete": True})
    with pytest.raises(ValueError, match="Refusing to replace"):
        write_json(path, {"complete": False})
    assert read_json(path) == {"complete": True}


def test_checkpoint_preserves_a_concurrent_writers_different_result(tmp_path, monkeypatch):
    path = tmp_path / "result.json"
    link = inventory.os.link

    def competing_writer(source, destination):
        Path(destination).write_bytes(b'{"winner":true}\n')
        link(source, destination)

    monkeypatch.setattr(inventory.os, "link", competing_writer)
    with pytest.raises(ValueError, match="Refusing to replace"):
        write_json(path, {"winner": False})
    assert read_json(path) == {"winner": True}


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
    result["limitations"].append("input_acquisition_truncated")
    candidate, = [c for c in summarize(result)["candidates"]
                   if c["amino_acids"] == "MPSRRRGGSSLIPGRMPILRFSVTTRYFSY"]
    assert "reconstruction_limit:input_acquisition_truncated" in candidate["uncertainty_flags"]
    assert candidate["rna_support"]["fragments"] == 15


@pytest.fixture
def empty_audit(tmp_path):
    """Real indexed BAM and pinned inputs for the offline acquisition/run/report path."""
    source = dict(id="source", url="original.bam", claims=[], selection=dict(cohort="tumor_candidate"))
    group = dict(breakends=[dict(contig="chr1", position=100, retained_side="left"),
                           dict(contig="chr1", position=200, retained_side="right")], nominations=["target"])
    manifest = dict(nominations={"target": dict(geometry_id="g")}, geometries={"g": group},
                    sources={"source": source})
    write_json(tmp_path / "inventory.json.gz", manifest)
    write_json(tmp_path / "inventory-pin.json", dict(sha256=digest(tmp_path / "inventory.json.gz")))
    reference = dict(models={}, assignments={"g": dict(transcripts=[], excluded_transcripts=[])},
                     inventory_sha256=digest(tmp_path / "inventory.json.gz"))
    write_json(tmp_path / "references.json.gz", reference)
    write_json(tmp_path / "references-pin.json", dict(sha256=digest(tmp_path / "references.json.gz")))
    bam = tmp_path / "empty.bam"
    header = dict(HD=dict(SO="coordinate"), SQ=[dict(SN="chr1", LN=1000)])
    with pysam.AlignmentFile(str(bam), "wb", header=header):
        pass
    pysam.index(str(bam))
    write_json(tmp_path / "sources" / "source" / "header.json",
               dict(source_id="source", source_identity=identity(source), status="ready", header=header))
    batch = dict(id="batch", geometries=["g"], regions=[["chr1", 0, 300]], window=100)
    write_json(tmp_path / "batches.json", [batch])
    return source, batch, bam


@pytest.mark.parametrize("status,limits", [("bounded", []), ("truncated", ["max_partner_queries"]),
                                         ("incomplete", ["max_rounds"])])
def test_bounded_inputs_remain_explicit_through_offline_report(tmp_path, monkeypatch, empty_audit, status, limits):
    source, batch, bam = empty_audit
    monkeypatch.setattr(acquisition, "source_file", lambda source: source)
    monkeypatch.setattr(acquisition, "extract_reads", lambda *args, **kwargs: fake_subset(bam, dict(status=status, limits=limits)))
    result = acquisition.acquire_batch(tmp_path, source, batch, None, osteosarc.RecoveryPolicy())
    expected = "acquired" if status == "bounded" else status
    assert result["status"] == expected
    reconstruct.run_batch(tmp_path, "source", "batch", reconstruct.PARAMETERS, 30, "engine")
    for view in ("forward", "reverse"):
        result = read_json(reconstruct.result_path(tmp_path, "source", "g", view))
        assert result["status"] == "no_candidate_paths"
        assert result["acquisition_status"] == expected
        assert result["input_limitations"] == acquisition.input_limitations(dict(status=expected, receipt=dict(limits=limits)))
        assert set(result["input_limitations"]) <= set(result["limitations"])
    outcome = report.report(tmp_path, tmp_path / "report", require_complete=True)
    assert outcome["all_pairs_accounted"]
    ledger = read_json(tmp_path / "report" / "coverage.json.gz")
    assert {row["acquisition_status"] for row in ledger["outcomes"]} == {expected}
    # This confirms jobs ran, not that incomplete inputs establish biological absence.
    assert ledger["counts"] == {"accounted": 1}


def test_absent_contig_keeps_other_targets_and_cannot_become_a_negative(tmp_path, monkeypatch, empty_audit):
    source, batch, bam = empty_audit
    batch["regions"].append(["chrMissing", 0, 100])
    captured = []
    monkeypatch.setattr(acquisition, "source_file", lambda source: source)

    def extract(source, regions, **kwargs):
        captured.extend(regions)
        return fake_subset(bam, dict(status="bounded"))

    monkeypatch.setattr(acquisition, "extract_reads", extract)
    result = acquisition.acquire_batch(tmp_path, source, batch, None, osteosarc.RecoveryPolicy())
    assert result["status"] == "acquired"
    assert result["unavailable_contigs"] == ["chrMissing"]
    assert [r.contig for r in captured] == ["chr1"]
    reconstruct.run_batch(tmp_path, "source", "batch", reconstruct.PARAMETERS, 30, "engine")
    result = read_json(reconstruct.result_path(tmp_path, "source", "g", "forward"))
    assert result["status"] == "no_candidate_paths"  # The unaffected target was still queried.
    assert "input_acquisition_missing_contigs" in result["limitations"]
    # The same missing contig at the nominated event is explicitly unassessable.
    manifest = load_inventory(tmp_path)
    manifest["geometries"]["g"]["breakends"][0]["contig"] = "chrMissing"
    monkeypatch.setattr(reconstruct, "load_inventory", lambda directory: manifest)
    for view in ("forward", "reverse"):
        reconstruct.result_path(tmp_path, "source", "g", view).unlink()
    reconstruct.run_batch(tmp_path, "source", "batch", reconstruct.PARAMETERS, 30, "engine")
    result = read_json(reconstruct.result_path(tmp_path, "source", "g", "forward"))
    assert result["status"] == "not_assessable"
    assert result["reason"] == "missing_or_ambiguous_contig"


def test_report_rejects_mixed_engines_and_corrupted_details(tmp_path, monkeypatch, empty_audit):
    source, batch, bam = empty_audit
    monkeypatch.setattr(acquisition, "source_file", lambda source: source)
    monkeypatch.setattr(acquisition, "extract_reads", lambda *args, **kwargs: fake_subset(bam, dict(status="bounded")))
    acquisition.acquire_batch(tmp_path, source, batch, None, osteosarc.RecoveryPolicy())
    reconstruct.run_batch(tmp_path, "source", "batch", reconstruct.PARAMETERS, 30, "engine")
    path = reconstruct.result_path(tmp_path, "source", "g", "reverse")
    original = read_json(path)
    changed = deepcopy(original)
    changed["request"]["engine_id"] = "different-engine"
    path.unlink()
    write_json(path, changed)
    with pytest.raises(ValueError, match="mix engine"):
        report.report(tmp_path, tmp_path / "report")
    path.unlink()
    write_json(path, original)
    (tmp_path / original["reconstruction"]["path"]).write_bytes(b"corrupt")
    with pytest.raises(ValueError, match="checksum mismatch"):
        report.report(tmp_path, tmp_path / "report")


@pytest.mark.parametrize("timeout", [False, True])
def test_deadline_excludes_persistence_and_timeout_has_no_partial_candidates(
        tmp_path, monkeypatch, empty_audit, timeout):
    source, batch, bam = empty_audit
    monkeypatch.setattr(acquisition, "source_file", lambda source: source)
    monkeypatch.setattr(acquisition, "extract_reads", lambda *args, **kwargs: fake_subset(bam, dict(status="bounded")))
    acquisition.acquire_batch(tmp_path, source, batch, None, osteosarc.RecoveryPolicy())
    armed = []
    monkeypatch.setattr(reconstruct.signal, "alarm", lambda seconds: armed.append(seconds))
    persist = reconstruct.write_json

    def check_persistence(path, value):
        assert armed[-1] == 0
        persist(path, value)

    monkeypatch.setattr(reconstruct, "write_json", check_persistence)
    if timeout:
        def interrupted_summary(result):
            assert armed[-1] == 30
            raise TimeoutError("interrupted summary")
        monkeypatch.setattr(reconstruct, "summarize", interrupted_summary)
    reconstruct.run_batch(tmp_path, "source", "batch", reconstruct.PARAMETERS, 30, "engine")
    for orientation in ("forward", "reverse"):
        result = read_json(reconstruct.result_path(tmp_path, "source", "g", orientation))
        if timeout:
            assert result["status"] == "reconstruction_timeout"
            assert "candidates" not in result and "reconstruction" not in result
            assert not (tmp_path / "reconstructions").exists()
        else:
            assert result["status"] == "no_candidate_paths"
            assert "reconstruction" in result


def test_insufficient_cache_space_leaves_work_pending(tmp_path, monkeypatch, empty_audit):
    _, batch, _ = empty_audit
    monkeypatch.setattr(acquisition, "planned_batches", lambda *args: [batch])
    monkeypatch.setattr(acquisition.shutil, "disk_usage", lambda path: SimpleNamespace(free=0))
    with pytest.raises(OSError, match="Acquisition paused"):
        acquisition.acquire(tmp_path, tmp_path / "cache", workers=1)
    assert not (tmp_path / "sources" / "source" / "batch.json").exists()
    with pytest.raises(ValueError, match="Audit is incomplete"):
        report.report(tmp_path, tmp_path / "report", require_complete=True)
    ledger = read_json(tmp_path / "report" / "coverage.json.gz")
    assert ledger["counts"] == {"pending_reconstruction": 1}


def test_acquisition_pins_shared_provenance_and_rejects_modified_assets(tmp_path, monkeypatch, empty_audit):
    from osteosarc.read_receipts import write_read_receipt

    source, batch, bam = empty_audit
    workspace = tmp_path / "cache" / "osteosarc"
    derived = workspace / "derived" / "subset"
    derived.mkdir(parents=True)
    for suffix in ("", ".bai"):
        shutil.copyfile(str(bam) + suffix, str(derived / "reads.bam") + suffix)
    receipt = dict(status="bounded", records=0, limits=[], acquisition=[
        dict(request=dict(filters=dict(query_names=["a", "b"]))) for _ in range(5)])
    write_read_receipt(derived / "receipt.json", receipt, workspace)
    subset = osteosarc.ReadSubset(derived / "reads.bam", derived / "reads.bam.bai", receipt)
    monkeypatch.setattr(acquisition, "source_file", lambda source: source)
    monkeypatch.setattr(acquisition, "extract_reads", lambda *args, **kwargs: subset)
    saved = acquisition.acquire_batch(tmp_path, source, batch, None, osteosarc.RecoveryPolicy())
    assert saved["receipt"] == dict(status="bounded", records=0, limits=[])
    assert len(saved["upstream_receipt_files"]) == 2
    with pytest.raises(ValueError, match="request drift"):
        acquisition.acquire_batch(tmp_path, source, batch, None, osteosarc.RecoveryPolicy(), timeout=301)
    reconstruct.run_batch(tmp_path, "source", "batch", reconstruct.PARAMETERS, 30, "engine")
    asset, = (workspace / "query-names").glob("*.txt")
    asset.write_text("corrupt\n")
    with pytest.raises(ValueError, match="provenance checksum"):
        acquisition.acquire_batch(tmp_path, source, batch, None, osteosarc.RecoveryPolicy())
    with pytest.raises(ValueError, match="provenance checksum"):
        reconstruct.run_batch(tmp_path, "source", "batch", reconstruct.PARAMETERS, 30, "engine")
    with pytest.raises(ValueError, match="provenance checksum"):
        report.report(tmp_path, tmp_path / "report")


def test_oversized_batch_is_partitioned_without_losing_targets_or_context(tmp_path, monkeypatch, empty_audit):
    source, _, bam = empty_audit
    manifest = load_inventory(tmp_path)
    manifest["geometries"]["other"] = deepcopy(manifest["geometries"]["g"])
    for end in manifest["geometries"]["other"]["breakends"]:
        end["position"] += 500
    manifest["nominations"]["other"] = dict(geometry_id="other")
    reference = read_json(tmp_path / "references.json.gz")
    reference["assignments"]["other"] = dict(transcripts=[], excluded_transcripts=[])
    for name in ("inventory.json.gz", "inventory-pin.json", "references.json.gz", "references-pin.json", "batches.json"):
        (tmp_path / name).unlink()
    write_json(tmp_path / "inventory.json.gz", manifest)
    write_json(tmp_path / "inventory-pin.json", dict(sha256=digest(tmp_path / "inventory.json.gz")))
    reference["inventory_sha256"] = digest(tmp_path / "inventory.json.gz")
    write_json(tmp_path / "references.json.gz", reference)
    write_json(tmp_path / "references-pin.json", dict(sha256=digest(tmp_path / "references.json.gz")))
    batch = acquisition.make_batch(manifest, reference, ["g", "other"], window=50)
    write_json(tmp_path / "batches.json", [batch])
    queried = []
    monkeypatch.setattr(acquisition, "source_file", lambda source: source)

    def extract(source, regions, **kwargs):
        queried.append(regions)
        if len(regions) > 1:
            raise acquisition.RecordLimitError("Acquisition exceeds record limit")
        return fake_subset(bam, dict(status="bounded"))

    monkeypatch.setattr(acquisition, "extract_reads", extract)
    result = acquisition.acquire_batch(tmp_path, source, batch, None, osteosarc.RecoveryPolicy(),
                                       manifest=manifest, reference=reference)
    assert result["status"] == "partitioned"
    assert len(queried) == 3
    leaves, receipts = acquisition.acquisition_tree(tmp_path, "source", batch["id"])
    assert len(leaves) == 2 and len(receipts) == 3
    outcome = reconstruct.run_batch(tmp_path, "source", batch["id"], reconstruct.PARAMETERS, 30, "engine")
    assert outcome["outcomes"] == {"no_candidate_paths": 4}
    assert report.report(tmp_path, tmp_path / "report", require_complete=True)["expected_pairs"] == 2
    # A successful sibling does not excuse a missing child receipt.
    child_id = next(iter(leaves))
    (tmp_path / "sources" / "source" / (child_id + ".json")).unlink()
    with pytest.raises(FileNotFoundError):
        report.report(tmp_path, tmp_path / "missing-child")
