"""Shared checks for packaged full-depth performance evidence."""

import json

from tests.data.osteosarc.expansion.inventory import digest


# Inputs and settings both sides of a collection pair must share. The packaging
# validator, collection_report.collect_runs, enforces the same keys.
COLLECTION_PAIR_IDENTITY_KEYS = (
    "source_bam_sha256", "receipt_sha256", "input_bam_sha256", "input_index_sha256", "inventory_sha256",
    "benchmark_sha256", "variant_id", "mode", "merge_overlapping_fragments", "allocation_instrumented")
PIPELINE_PAIR_IDENTITY_KEYS = (
    "source_bam_sha256", "receipt_sha256", "input_bam_sha256", "input_index_sha256", "inventory_sha256",
    "reference_manifest_sha256", "benchmark_sha256", "variant_id", "mode")
COLLECTION_RESULT_KEYS = ("locus_reads_sha256", "allele_reads_sha256", "locus_read_count", "allele_read_count", "counts")
PIPELINE_RESULT_KEYS = ("rna_sha256", "proteins_sha256", "counts", "public_protein_count", "checked_protein_count")


def load_pinned_package(results, generator):
    """Verify a package's manifest and generator digest, then load its reports."""
    manifest = json.loads((results / "manifest.json").read_text())
    assert set(manifest) == {"files", "generator_sha256"}
    for name, checksum in manifest["files"].items():
        assert digest(results / name) == checksum
    # The shipped validator must be the one that reproduces this package; a
    # validator change requires repackaging or re-pinning in the same commit.
    assert manifest["generator_sha256"] == digest(generator)
    return {name: json.loads((results / (name + ".json")).read_text()) for name in ("collection", "pipeline")}


def changed_sources(before, after):
    """Source files whose digest differs between two runs, including added or removed files."""
    left, right = before["identity"]["source_files"], after["identity"]["source_files"]
    return {path for path in set(left) | set(right) if left.get(path) != right.get(path)}
