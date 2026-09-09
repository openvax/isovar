"""Completeness, missing-state, attribution and acquisition-cache safeguards."""

from copy import deepcopy
import json

import pytest

from tests.data.osteosarc.expansion.acquire import survey_header
from tests.data.osteosarc.expansion.discover import list_prefix
from tests.data.osteosarc.expansion.inventory import digest, fetch_snapshot
from tests.data.osteosarc.expansion.report import (
    compact_row, missing_row, protein_supporting_read_names, source_disposition, summarize,
)


@pytest.mark.parametrize("header,indexes,receipts,expected", [
    (None, [], [], "header_unavailable"),
    ({"status": "error"}, [], [], "header_unavailable"),
    ({"status": "ok", "coordinate_status": "unmapped_reference"}, [], [], "no_genomic_coordinates"),
    ({"status": "ok", "assembly": None}, [], [], "unresolved_reference"),
    ({"status": "ok", "assembly": "GRCh38"}, [], [], "missing_index"),
    ({"status": "ok", "assembly": "GRCh38"}, ["index"], [], "acquisition_not_attempted"),
    ({"status": "ok", "assembly": "GRCh38"}, ["index"], [{"status": "error"}], "acquisition_failed"),
    ({"status": "ok", "assembly": "GRCh38"}, ["index"], [{"status": "ok"}], "audit_not_completed"),
    ({"status": "ok", "assembly": "GRCh38"}, ["index"],
     [{"status": "ok", "variant_ids": ["nuclear"]}, {"status": "error", "variant_ids": ["MT"]}], "acquisition_failed"),
])
def test_missing_matrix_cell_is_never_biological_zero(header, indexes, receipts, expected):
    result = missing_row(dict(source_id="source", listed_indexes=indexes), dict(variant_id="MT"), header, receipts)
    assert result["status"] == expected
    for mode in ("defaults", "primary_only"):
        assert result[mode]["counts"] is None
        assert result[mode]["status"] == "unavailable"
    assert result["independent_primary"] is None


def test_report_keeps_support_names_all_checks_and_rank_order():
    original = dict(defaults=dict(proteins=[dict(amino_acids="SHORT", supporting_read_names=["b", "a"], checks=[2, 1]),
                                           dict(amino_acids="LONGER", supporting_read_names=["a"], checks=[3])]))
    untouched = deepcopy(original)
    row = compact_row(original)
    assert original == untouched
    proteins = row["defaults"]["proteins"]
    assert [p["amino_acids"] for p in proteins] == ["SHORT", "LONGER"]
    assert protein_supporting_read_names(row, proteins[0]) == ["a", "b"]
    assert protein_supporting_read_names(row, proteins[1]) == ["a"]
    assert proteins[0]["checks"] == [1, 2]
    assert len(proteins[0]["supporting_read_names_sha256"]) == 64


@pytest.mark.parametrize("indices", [[-1], [2], [1, 0], [0, 0], [False], ["0"]])
def test_report_rejects_invalid_support_name_indices(indices):
    row = compact_row(dict(defaults=dict(proteins=[dict(supporting_read_names=["a", "b"], checks=[])])))
    protein = row["defaults"]["proteins"][0]
    protein["supporting_read_name_indices"] = indices
    with pytest.raises(ValueError):
        protein_supporting_read_names(row, protein)


def test_report_rejects_tampered_support_name_table():
    row = compact_row(dict(defaults=dict(proteins=[dict(supporting_read_names=["original"], checks=[])])))
    row["supporting_read_name_table"] = ["changed"]
    with pytest.raises(ValueError, match="checksum mismatch"):
        protein_supporting_read_names(row, row["defaults"]["proteins"][0])


def test_source_classification_requires_evidence_not_reference_name_alone():
    source = dict(key="unknown.bam", scope="unresolved_assay")
    assert source_disposition(source, None, {})[0] == "unresolved_assay"
    assert source_disposition(source, None, {"unknown.bam": "WGS"})[0] == "excluded_dna_alignment"
    hla = dict(status="ok", header=dict(SQ=[dict(SN="HLA00001", LN=546)]))
    assert source_disposition(dict(source, scope="rna_candidate"), hla, {})[0] == "excluded_receptor_or_targeted_alignment"
    genomic = dict(status="ok", header=dict(SQ=[dict(SN="chr1", LN=248956422)]))
    assert source_disposition(source, genomic, {})[0] == "unresolved_assay"
    dragen = dict(status="ok", header=dict(PG=[dict(ID=" DRAGEN SW build", CL="Mapper.rna-mode: 0 ")]))
    assert source_disposition(dict(source, key="kamil/basespace/results/source.bam"), dragen, {})[0] == "excluded_dna_alignment"


def test_summary_retains_timeout_counts_without_claiming_protein():
    result = dict(status="resource_limit", outcome="audit_timeout", counts=dict(reads=dict(ref=1, alt=100, other=2)))
    rows = [dict(variant_id="variant", status="audited", defaults=result, primary_only=result)]
    summary, = summarize(rows, [dict(variant_id="variant", gene="GENE")])
    assert summary["defaults"]["products_with_alt"] == 1
    assert summary["defaults"]["products_with_ranked_protein"] == 0
    assert summary["defaults"]["outcomes"] == {"audit_timeout": 1}


def test_header_cache_refuses_source_or_checksum_drift(tmp_path, monkeypatch):
    directory = tmp_path / "alignments/source"
    directory.mkdir(parents=True)
    header = directory / "header.sam"
    header.write_text("@HD\tVN:1.6\n")
    receipt = dict(source_url="https://example/original.bam", status="ok", header_sha256=digest(header))
    (directory / "header.json").write_text(json.dumps(receipt))
    monkeypatch.setattr("subprocess.run", lambda *a, **k: pytest.fail("Cached header accessed network"))
    source = dict(source_id="source", url=receipt["source_url"])
    assert survey_header(source, tmp_path) == receipt
    with pytest.raises(ValueError, match="source mismatch"):
        survey_header(dict(source, url="https://different/"), tmp_path)
    header.write_text("corrupt")
    with pytest.raises(ValueError, match="checksum mismatch"):
        survey_header(source, tmp_path)


@pytest.mark.parametrize("state", ["receipt_only", "partial_only", "snapshot_only"])
def test_incomplete_snapshot_fails_closed(tmp_path, monkeypatch, state):
    path = tmp_path / "snapshot"
    if state == "receipt_only":
        path.with_name("snapshot.receipt.json").write_text("{}")
    elif state == "partial_only":
        path.with_name("snapshot.partial").write_text("incomplete")
    else:
        path.write_text("unreceipted")
    monkeypatch.setattr("subprocess.run", lambda *a, **k: pytest.fail("Invalid cache accessed network"))
    with pytest.raises((ValueError, FileNotFoundError)):
        fetch_snapshot("https://example/", path)


@pytest.mark.parametrize("mode", ["complete", "repeat_token", "wrong_prefix", "missing_token"])
def test_paginated_discovery_is_bounded_and_fail_closed(tmp_path, monkeypatch, mode):
    from tests.data.osteosarc.expansion import discover
    calls = []

    def fetch(url, path):
        calls.append(url)
        prefix = "other/" if mode == "wrong_prefix" else "rna/"
        truncated = len(calls) == 1 or mode != "complete"
        token = "" if mode == "missing_token" else "<NextContinuationToken>same</NextContinuationToken>"
        path.write_text('<ListBucketResult xmlns="http://s3.amazonaws.com/doc/2006-03-01/">'
                        f"<Prefix>{prefix}</Prefix><IsTruncated>{str(truncated).lower()}</IsTruncated>{token}"
                        "</ListBucketResult>")
        return dict(url=url, sha256=digest(path))

    monkeypatch.setattr(discover, "fetch_snapshot", fetch)
    if mode == "complete":
        records, receipts = list_prefix("rna/", tmp_path)
        assert records == [] and len(receipts) == 2
        assert "continuation-token=same" in calls[1]
    else:
        with pytest.raises(ValueError):
            list_prefix("rna/", tmp_path)
        assert len(calls) <= 2
