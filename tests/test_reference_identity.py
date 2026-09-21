"""Offline identities follow verified reference contents, not paths or labels."""

import gzip
from hashlib import sha256
import json
from pathlib import Path
import shutil

import pytest
from varcode import Variant

from tests.data.osteosarc.expansion.references import reference_genome as expansion_genome
from tests.osteosarc_protein_helpers import REFERENCE, reference_genome as protein_genome
from tests.reference_identity import REFERENCE_FILES, reference_dataset_identity


STRESS = REFERENCE.parent / "expansion/stress-corpus/reference"


def write_assets(directory, contents, mtime=0):
    directory.mkdir(exist_ok=True)
    checksums = {}
    for name, data in contents.items():
        archive = gzip.compress(data, mtime=mtime)
        (directory / name).write_bytes(archive)
        checksums[name] = sha256(archive).hexdigest()
    return checksums


def test_identity_survives_relocation_and_recompression(tmp_path):
    contents = dict(zip(REFERENCE_FILES, (b"annotation", b"transcripts", b"proteins")))
    original, moved = tmp_path / "original", tmp_path / "moved"
    pins = write_assets(original, contents)
    identity = reference_dataset_identity(original, pins)
    shutil.copytree(original, moved)
    assert reference_dataset_identity(moved, pins) == identity
    repacked = write_assets(moved, contents, mtime=1234)
    assert repacked != pins
    assert reference_dataset_identity(moved, repacked) == identity


@pytest.mark.parametrize("name", REFERENCE_FILES)
def test_each_reference_asset_changes_identity_and_requires_updated_pin(tmp_path, name):
    contents = dict(zip(REFERENCE_FILES, (b"annotation", b"transcripts", b"proteins")))
    pins = write_assets(tmp_path, contents)
    identity = reference_dataset_identity(tmp_path, pins)
    contents[name] += b" changed"
    updated = write_assets(tmp_path, contents)
    with pytest.raises(ValueError, match="checksum mismatch"):
        reference_dataset_identity(tmp_path, pins)
    assert reference_dataset_identity(tmp_path, updated) != identity


def test_reference_roles_and_required_pins_are_not_interchangeable(tmp_path):
    contents = dict(zip(REFERENCE_FILES, (b"annotation", b"transcripts", b"proteins")))
    pins = write_assets(tmp_path, contents)
    identity = reference_dataset_identity(tmp_path, pins)
    contents[REFERENCE_FILES[1]], contents[REFERENCE_FILES[2]] = contents[REFERENCE_FILES[2]], contents[REFERENCE_FILES[1]]
    pins = write_assets(tmp_path, contents)
    assert reference_dataset_identity(tmp_path, pins) != identity
    for name in REFERENCE_FILES:
        with pytest.raises(ValueError, match="Missing reference asset checksum"):
            reference_dataset_identity(tmp_path, {k: v for k, v in pins.items() if k != name})


@pytest.mark.parametrize("kind", ["protein", "expansion"])
def test_offline_genome_relocation_and_metadata_do_not_change_identity(tmp_path, kind):
    source = REFERENCE if kind == "protein" else STRESS
    moved = tmp_path / "moved"
    shutil.copytree(source, moved)
    manifest_path = moved / ("protein_reference_manifest.json" if kind == "protein" else "manifest.json")
    manifest = json.loads(manifest_path.read_text())
    manifest["dataset_identity"] = "a different descriptive label"
    manifest["note"] = "moved to another directory"
    manifest_path.write_text(json.dumps(manifest, indent=4))
    make = (lambda directory, cache: protein_genome(cache, directory)) if kind == "protein" else expansion_genome
    first = make(source, tmp_path / "cache")
    second = make(moved, tmp_path / "cache")
    assert first.reference_name == second.reference_name
    assert first.contigs() == second.contigs()
    assert first.transcript_ids() == second.transcript_ids()


@pytest.mark.parametrize("kind", ["protein", "expansion"])
@pytest.mark.parametrize("name", REFERENCE_FILES)
def test_tampering_is_rejected_before_indexing(tmp_path, monkeypatch, kind, name):
    source = REFERENCE if kind == "protein" else STRESS
    moved = tmp_path / "moved"
    shutil.copytree(source, moved)
    (moved / name).write_bytes(gzip.compress(b"changed reference", mtime=0))

    def forbidden_index(*args, **kwargs):
        pytest.fail("An unverified reference reached indexing")

    monkeypatch.setattr("pyensembl.Genome.index", forbidden_index)
    with pytest.raises(ValueError, match="checksum mismatch"):
        if kind == "protein":
            protein_genome(tmp_path / "cache", moved)
        else:
            expansion_genome(moved, tmp_path / "cache")


@pytest.mark.parametrize("order", [("protein", "stress"), ("stress", "protein")])
def test_six_transcript_and_stress_annotations_are_independent_in_both_orders(tmp_path, monkeypatch, order):
    if hasattr(Variant, "_reference_name_to_valid_contig_names"):
        monkeypatch.setattr(Variant, "_reference_name_to_valid_contig_names", {})
    genomes = {}
    for name in order:
        genome = (protein_genome(tmp_path / "cache") if name == "protein"
                  else expansion_genome(STRESS, tmp_path / "cache"))
        genomes[name] = genome
        if name == "protein":
            # This subset has chromosome 1 but not the stress dataset's chr5.
            assert "PIP5K1A" in Variant("1", 151242178, "AG", "A", ensembl=genome).gene_names
        else:
            assert "ACSL6" in Variant("5", 131988563, "G", "GTA", ensembl=genome).gene_names
    assert genomes["protein"].reference_name != genomes["stress"].reference_name
