"""Only test-selected, byte-identical original Sid records reach the package."""

from collections import Counter
from copy import deepcopy
from importlib.metadata import version
from pathlib import Path
from types import SimpleNamespace
import subprocess
import sys

import pysam
import pytest

from isovar import sid_data
from tests.data.osteosarc import bundle as recipe_compiler
from tests.data.osteosarc.bundle import fixture_inputs


@pytest.fixture(scope="module")
def packaged():
    return sid_data.read_json(sid_data.RECIPE), sid_data.read_json(sid_data.BUNDLE)


def test_recipe_covers_exactly_the_consumed_original_records(packaged):
    recipe, bundle = packaged
    fixtures = list(fixture_inputs())
    assert set(recipe["fixtures"]) == {f["name"] for f in fixtures}
    for fixture in fixtures:
        selected = recipe["fixtures"][fixture["name"]]
        source = bundle["sources"][selected["source"]]
        assert source["asset"]["url"] == fixture["url"]
        assert selected["tests"] == fixture["tests"]
        assert [source["records"][key] for key in selected["records"]] == fixture["records"]
    assert sid_data.verify(recipe, bundle)["unique_records"] < sum(len(f["records"]) for f in fixtures)


def test_every_bundled_record_is_required_and_unaltered(packaged):
    recipe, bundle = packaged
    result = sid_data.verify(recipe, bundle)
    assert result["sources"] == 39
    assert result["fixtures"] == 311
    assert all(source["regions"] for source in recipe["sources"].values())
    assert {s["assembly"] for s in recipe["sources"].values()} == {"GRCh37", "GRCh38"}
    assert {r["reference_length"] for s in recipe["sources"].values() for r in s["regions"]
            if r["contig"] in ("M", "MT", "chrM", "chrMT")} == {16569, 16571}
    for source in bundle["sources"].values():
        header = pysam.AlignmentHeader.from_dict(source["header"])
        for line in source["records"].values():
            assert pysam.AlignedSegment.fromstring(line, header).to_string() == line


def test_offline_export_preserves_record_order_multiplicity_and_missing_qual(packaged, tmp_path, monkeypatch):
    recipe, bundle = packaged
    monkeypatch.setattr(subprocess, "run", lambda *a, **kw: pytest.fail("Unexpected subprocess"))
    monkeypatch.setattr(sid_data, "open_dataset", lambda *a, **kw: pytest.fail("Unexpected acquisition"))
    # Include mixed-platform clipping examples, hard-clipped split paths,
    # native mitochondrial references and a duplicate-containing legacy fixture.
    names = [n for n in recipe["fixtures"] if n.startswith(("read_ends/", "chimeric/"))]
    names += ["osteosarc/bulk_star_t0.sam.gz"]
    names += [n for n in recipe["fixtures"] if "MT_ND5" in n and n.endswith(".bam")]
    for index, name in enumerate(names):
        path = tmp_path / (str(index) + ".sam")
        sid_data.export_fixture(name, path, recipe=recipe, bundle=bundle)
        fixture = recipe["fixtures"][name]
        with pysam.AlignmentFile(path) as handle:
            assert [sid_data.sam_digest(r.to_string()) for r in handle] == fixture["records"]
        with pytest.raises(FileExistsError):
            sid_data.export_fixture(name, path, recipe=recipe, bundle=bundle)


def test_unneeded_or_changed_records_are_rejected(packaged):
    recipe, original = packaged
    bundle = deepcopy(original)
    source = next(iter(bundle["sources"].values()))
    key = next(iter(source["records"]))
    source["records"][key] += "\tZZ:Z:changed"
    with pytest.raises(ValueError, match="checksum"):
        sid_data.verify(recipe, bundle)
    source["records"].pop(key)
    with pytest.raises(ValueError, match="differ"):
        sid_data.verify(recipe, bundle)


def test_selector_requires_true_duplicate_multiplicity_and_preserves_unknown_quality():
    header = pysam.AlignmentHeader.from_dict({"SQ": [{"SN": "chr1", "LN": 1000}]})
    line = "r\t0\tchr1\t10\t60\t4M\t*\t0\t0\tACGT\t*\tZZ:Z:original"
    record = pysam.AlignedSegment.fromstring(line, header)
    key = sid_data.sam_digest(line)
    with pytest.raises(ValueError, match="missing 1"):
        sid_data.select_records([record], Counter({key: 2}))
    assert sid_data.select_records([record, record], Counter({key: 2})) == {key: line}
    changed = pysam.AlignedSegment.fromstring(line.replace("ACGT", "ACGA"), header)
    with pytest.raises(ValueError, match="missing"):
        sid_data.select_records([changed, changed], Counter({key: 2}))


def test_repeated_observation_references_do_not_fabricate_source_duplicates():
    recipe = {"fixtures": {
        "raw": {"source": "source", "records": ["a", "a", "b"]},
        "other-use": {"source": "source", "records": ["a", "b"]},
        "cross-references": {"source": "source", "records": ["b"] * 10, "record_references": True},
    }}
    assert sid_data.required_counts(recipe, "source") == Counter(a=2, b=1)


def test_reduced_header_keeps_read_group_program_chain_and_sa_target():
    header = {"SQ": [{"SN": c, "LN": 1000} for c in ("chr1", "chr2", "chr3")],
              "RG": [{"ID": "library", "PG": "aligner"}, {"ID": "unused"}],
              "PG": [{"ID": "basecaller"}, {"ID": "aligner", "PP": "basecaller"}, {"ID": "unused"}]}
    line = "r\t0\tchr1\t10\t60\t4M\t*\t0\t0\tACGT\t*\tRG:Z:library\tSA:Z:chr2,20,+,4M,60,0;"
    reduced = sid_data.minimal_header(header, [line])
    assert [r["SN"] for r in reduced["SQ"]] == ["chr1", "chr2"]
    assert reduced["RG"] == [{"ID": "library", "PG": "aligner"}]
    assert reduced["PG"] == header["PG"][:2]


@pytest.fixture
def osteosarc():
    import osteosarc
    return osteosarc


def test_compiled_recipe_records_installed_acquisition_version(packaged, monkeypatch, osteosarc):
    recipe, bundle = packaged
    fixture = next(fixture_inputs())
    selected = recipe["fixtures"][fixture["name"]]
    original_source = recipe["sources"][selected["source"]]
    asset = SimpleNamespace(**original_source["asset"])
    dataset = SimpleNamespace(
        id=recipe["snapshot_id"], asset=lambda url: asset,
        inspect_alignment=lambda key: SimpleNamespace(
            assembly=original_source["assembly"], header=bundle["sources"][selected["source"]]["header"]))
    monkeypatch.setattr(recipe_compiler, "fixture_inputs", lambda: iter([fixture]))

    compiled = recipe_compiler.compile_recipe(dataset)

    assert compiled["osteosarc_version"] == version("osteosarc")
    assert compiled["snapshot_id"] == recipe["snapshot_id"]
    assert compiled["fixtures"] == {fixture["name"]: selected}
    assert compiled["sources"][selected["source"]]["asset"] == original_source["asset"]
    # Recompilation must not rewrite the historical packaged provenance.
    assert sid_data.read_json(sid_data.RECIPE) == recipe


def test_regeneration_delegates_explicit_intervals_then_drops_background(tmp_path, osteosarc):
    path = tmp_path / "regional.sam"
    line = "required\t0\tchr1\t10\t60\t4M\t*\t0\t0\tACGT\t*"
    path.write_text("@SQ\tSN:chr1\tLN:1000\n" + line + "\n" + line.replace("required", "background") + "\n")
    asset = SimpleNamespace(id="asset", key="key", url="https://example/key.bam", size=100, modified=1)
    key = sid_data.sam_digest(line)
    recipe = {"sources": {"source": {"asset": sid_data.asset_identity(asset),
              "regions": [{"contig": "chr1", "start": 9, "end": 10, "assembly": "GRCh38"}]}},
              "fixtures": {"fixture": {"source": "source", "records": [key]}}}
    calls = []

    def extract(actual_asset, regions):
        assert actual_asset is asset
        calls.append(regions)
        return SimpleNamespace(receipt={"records": 2}, open=lambda: pysam.AlignmentFile(path))

    dataset = SimpleNamespace(id="snapshot", asset=lambda key: asset, extract_reads=extract)
    destination = tmp_path / "selected"
    sid_data.generate(recipe, destination, dataset)
    assert calls == [[osteosarc.Region("chr1", 9, 10, "GRCh38")]]
    selected = sid_data.read_json(destination / "source.json.gz")
    assert selected["records"] == {key: line}
    sid_data.generate(recipe, destination, dataset)
    assert len(calls) == 1  # Reuse only a validated selection.
    asset.size = 101
    with pytest.raises(ValueError, match="asset changed"):
        sid_data.generate(recipe, destination, dataset)


def test_sam_coordinate_conversion_and_explicit_mitochondrial_identity(osteosarc):
    assert sid_data.sam_regions(["chrM:12995-12996"], "GRCh37", {"chrM": 16571}) == [
        osteosarc.Region("chrM", 12994, 12996, "GRCh37", reference_length=16571)]
    assert sid_data.sam_regions(["chr1:1-1"], "GRCh38") == [osteosarc.Region("chr1", 0, 1, "GRCh38")]
    with pytest.raises(ValueError):
        sid_data.sam_regions([], "GRCh38")
    with pytest.raises(ValueError):
        sid_data.sam_regions(["chr1"], "GRCh38")


def test_metadata_path_cannot_download_a_whole_alignment():
    with pytest.raises(ValueError, match="Whole alignments"):
        sid_data.fetch_metadata("https://example/source.bam?download=1")


def test_snapshot_bound_metadata_and_index_paths(tmp_path, monkeypatch):
    path = tmp_path / "source.bai"
    path.write_bytes(b"index bytes")
    asset = SimpleNamespace(id="index-id")
    url = "https://sid-sijbrandij-osteosarc-dataset.s3.us-west-2.amazonaws.com/source.bam.bai"
    downloaded = []
    dataset = SimpleNamespace(id="snapshot", manifest={"sources": {}}, asset=lambda u: asset,
                              download=lambda a: downloaded.append(a) or path)
    monkeypatch.setattr(sid_data, "open_dataset", lambda: dataset)
    actual_path, receipt = sid_data.fetch_metadata(url)
    assert actual_path == path
    assert receipt["asset_id"] == asset.id and receipt["snapshot_id"] == "snapshot"
    assert downloaded == [asset]
    dataset.manifest["sources"]["metadata"] = dict(url="https://example/metadata", sha256="pinned")
    dataset.cache = SimpleNamespace(path=lambda r: path)
    assert sid_data.fetch_metadata("https://example/metadata")[1]["sha256"] == "pinned"
    assert len(downloaded) == 1


def test_existing_expansion_api_uses_osteosarc_and_keeps_receipts(tmp_path, monkeypatch, osteosarc):
    from tests.data.osteosarc.expansion import acquire
    header = {"SQ": [{"SN": "chr1", "LN": 248956422}, {"SN": "chr2", "LN": 242193529}]}
    header_path = tmp_path / "original.sam"
    header_path.write_text(str(pysam.AlignmentHeader.from_dict(header)))
    calls = []
    dataset = SimpleNamespace(inspect_alignment=lambda *a, **kw: SimpleNamespace(
        header=header, path=header_path, receipt={"header": "receipt"}))
    monkeypatch.setattr(acquire, "open_dataset", lambda: dataset)

    def extract(url, regions, assembly, output, **kwargs):
        calls.append((url, regions, assembly, kwargs["reference_lengths"]))
        output.write_bytes(b"bounded BAM")
        Path(str(output) + ".bai").write_bytes(b"bounded index")
        return SimpleNamespace(receipt={"records": 1, "request": {"snapshot_id": "snapshot"}})

    monkeypatch.setattr(acquire, "extract_regions", extract)
    source = dict(source_id="id", url="https://source", listed_indexes=["source.bai"])
    variants = [dict(chrom="chr1", pos=10, ref="A", variant_id="v")]
    result = acquire.acquire_regions(source, variants, tmp_path, "GRCh38")
    assert result["status"] == "ok" and result["region_records"] == 1
    assert result["osteosarc"]["request"]["snapshot_id"] == "snapshot"
    assert calls == [(source["url"], ["chr1:9-11"], "GRCh38", {"chr1": 248956422, "chr2": 242193529})]
    assert acquire.acquire_regions(source, variants, tmp_path, "GRCh38") == result
    assert len(calls) == 1
