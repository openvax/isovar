"""Offline integrity and exhaustive count/state checks for the pinned full audit.

The large support-name/uncapped-protein artifact is checksum/CRC streamed;
routine CI does not load millions of support indices in every test worker.
Original-read fixtures separately exercise the real translation pipeline.
"""

from collections import Counter
import gzip
import json
from pathlib import Path

import pytest

from tests.data.osteosarc.expansion.inventory import digest


DIRECTORY = Path(__file__).parent / "data/osteosarc/expansion/audit"


@pytest.fixture(scope="session")
def full_audit_index():
    with gzip.open(DIRECTORY / "counts-index.json.gz", "rt") as handle:
        return json.load(handle)


def test_checked_in_full_audit_checksums_and_compressed_integrity():
    manifest = json.loads((DIRECTORY / "manifest.json").read_text())
    assert manifest["row_count"] == 7216 and manifest["incomplete_row_count"] == 0
    for name, checksum in manifest["files"].items():
        assert digest(DIRECTORY / name) == checksum
        if name.endswith(".gz"):
            with gzip.open(DIRECTORY / name, "rb") as handle:
                while handle.read(1024 * 1024):
                    pass  # Reading through EOF validates the gzip CRC/length.


def test_full_matrix_has_every_variant_for_every_eligible_product(full_audit_index):
    data = full_audit_index
    assert len(data["variant_ids"]) == len(set(data["variant_ids"])) == 44
    assert len(data["source_ids"]) == len(set(data["source_ids"])) == 164
    keys = [(r["source_id"], r["variant_id"]) for r in data["rows"]]
    assert len(keys) == len(set(keys)) == 7216
    assert set(keys) == {(s, v) for s in data["source_ids"] for v in data["variant_ids"]}
    sources = data["sources"]
    assert len(sources) == len({s["source_id"] for s in sources}) == 841
    assert Counter(s["scope"] for s in sources) == dict(
        rna_candidate=164, excluded_dna_alignment=108, excluded_receptor_or_targeted_alignment=558,
        excluded_transcript_coordinate_alignment=10, unresolved_assay=1)
    assert {s["source_id"] for s in sources if s["scope"] == "rna_candidate"} == set(data["source_ids"])
    assert all(not s["published_md5_verified_against_full_remote_file"] for s in sources)
    assert all("CL" not in p for s in sources for p in s.get("programs", []))


def test_full_matrix_unavailable_and_timeout_cells_are_not_zero_support(full_audit_index):
    rows = full_audit_index["rows"]
    for mode in ("defaults", "primary_only"):
        assert sum(r[mode].get("counts") is not None for r in rows) == 6202
        missing = [r[mode] for r in rows if r[mode].get("counts") is None]
        assert Counter(m["outcome"] for m in missing) == dict(acquisition_failed=442, missing_index=220,
                                                            no_genomic_coordinates=352)
        assert all(m["status"] == "unavailable" for m in missing)
        timeouts = [r for r in rows if r[mode]["outcome"] == "audit_timeout"]
        assert len(timeouts) == 5  # Pinned outcomes on this audit host/budget, not a timing assertion for CI.
        for row in timeouts:
            assert row["variant_id"] == "MT_ND5-chrM-12994"
            assert row[mode]["status"] == "resource_limit"
            assert row[mode]["counts"]["reads"]["alt"] > 0
            assert row[mode]["stages"]["stage"] in ("rna_sequence_creation", "independent_protein_validation")
            assert row[mode]["error"]["type"] == "AuditTimeout"
    stages = Counter(r[m]["stages"]["stage"] for r in rows for m in ("defaults", "primary_only")
                     if r[m]["outcome"] == "audit_timeout")
    assert stages == dict(rna_sequence_creation=9, independent_protein_validation=1)


def test_all_completed_count_units_and_ranked_checks_are_consistent(full_audit_index):
    rows = full_audit_index["rows"]
    for row in rows:
        for mode in ("defaults", "primary_only"):
            value = row[mode]
            assert value["returned_validation_errors"] == value["uncapped_validation_errors"] == 0
            assert value["returned_protein_count"] <= value["uncapped_protein_count"]
            counts = value.get("counts")
            if counts is None:
                continue
            for group in ("ref", "alt", "other"):
                assert counts["reads"][group] >= counts["read_objects"][group] >= counts["template_names"][group] >= 0
            assert max(counts["template_names"].values()) <= counts["all_template_names"] <= sum(counts["template_names"].values())
    assert sum(r["defaults"]["uncapped_protein_count"] for r in rows) == 42049
    assert sum(r["primary_only"]["uncapped_protein_count"] for r in rows) == 41864
    assert sum(s["defaults"]["products_with_alt"] > 0 for s in full_audit_index["summary"]) == 43
    assert sum(s["defaults"]["products_with_ranked_protein"] > 0 for s in full_audit_index["summary"]) == 39


def test_full_matrix_native_mitochondrial_coordinates_and_alias_recovery(full_audit_index):
    rows = {r["source_id"]: r for r in full_audit_index["rows"] if r["variant_id"] == "MT_ND5-chrM-12994"}
    bulk = rows["f30f618fb76a0e49"]
    assert bulk["defaults"]["counts"]["reads"] == dict(ref=10, alt=27, other=1)
    assert bulk["normalized_variant"]["contig"] == "MT" and bulk["native_variant"]["chrom"] == "chrM"
    assert rows["fbeb5d8a415e4d18"]["native_variant"]["pos"] == 12995
    assert rows["238f917202dd9997"]["native_variant"]["pos"] == 12994
    assert rows["fbeb5d8a415e4d18"]["native_variant"]["assembly"] == "GRCh37"
    inventories = full_audit_index["native_variant_inventories"]
    assert len(inventories) == 3
    assert all(len(v["variants"]) == 44 for v in inventories.values())
    assert all(len(v["sha256"]) == 64 for v in inventories.values())
    for run in full_audit_index["audit_runs"].values():
        assert run["software"]["isovar"] == "1.8.1"
        assert run["time_limit_seconds_per_mode"] == 120
        assert run["settings"]["protein_sequence_creator"]["max_protein_sequences_per_variant"] == 1
