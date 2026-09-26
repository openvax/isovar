"""Isovar's Sid test reads come from openvax-v1, record for record."""

import json
from pathlib import Path
from types import SimpleNamespace

import pysam
import pytest

from isovar import sid_data

DATA = Path(__file__).parent / "data"


@pytest.fixture(scope="module")
def members():
    return sid_data.members()


def _records(path):
    with pysam.AlignmentFile(str(path), check_sq=False) as handle:
        return handle.header, [r.to_string() for r in handle]


def test_every_read_file_is_exported_in_its_original_format(members):
    files = [name for name in members if "#" not in name]
    assert len(files) == 124
    for name in files:
        path = sid_data.path(name)
        # osteosarc names each export <member>.<format>, keeping the original suffix.
        assert path.name == "%s.%s" % (Path(name).name, sid_data._format(name))
        if name.endswith(".bam"):
            assert Path(str(path) + ".bai").exists()
        _records(path)
    manifest = json.loads((DATA / "osteosarc/manifest.json").read_text())
    for data in manifest["datasets"].values():
        _, records = _records(sid_data.path("osteosarc/" + data["file"]))
        assert len(records) == data["fixture_records"]


def test_embedded_and_selected_reads_match_the_bundle(members, tmp_path):
    from osteosarc import check_fixtures
    selections = [name for name in members if "#" in name]
    assert len(selections) == 221
    fixtures = {}
    for name in selections:
        path, _, selector = name.partition("#")
        if selector.startswith("/"):
            fixtures[sid_data.PREFIX + name] = {"json": str(DATA / path), "pointer": selector}
        else:
            # One read of a SAM file; ":" in its name is "-" in the member name.
            with open(DATA / path) as handle:
                lines = [line for line in handle if line.startswith("@") or
                         line.split("\t")[0] in (selector, selector.replace("-", ":"))]
            single = tmp_path / ("%d.sam" % len(fixtures))
            single.write_text("".join(lines))
            fixtures[sid_data.PREFIX + name] = str(single)
    assert check_fixtures(sid_data.bundle(), fixtures) == {}


@pytest.mark.parametrize("name", ["osteosarc/no-such-file.bam", "osteosarc", "", "/etc/hosts",
                                  "../README.md", "fusions/corpus/TPST1--CRCP-T1.input.json.gz#/original_records"])
def test_path_accepts_only_isovar_read_files(name):
    with pytest.raises(KeyError):
        sid_data.path(name)


def test_export_never_overwrites_and_needs_a_known_format(tmp_path):
    output = tmp_path / "chimeric.bam"
    sid_data.export("chimeric/osteosarc-ont.sam", output)
    assert len(_records(output)[1]) == 2 and Path(str(output) + ".bai").exists()
    with pytest.raises(FileExistsError):
        sid_data.export("chimeric/osteosarc-ont.sam", output)
    other = tmp_path / "other.bam"
    Path(str(other) + ".bai").write_text("existing index")
    with pytest.raises(FileExistsError):
        sid_data.export("chimeric/osteosarc-ont.sam", other)
    for name in ("reads.SAM", "reads.cram", "reads.txt"):
        with pytest.raises(ValueError, match="bam, .sam or .sam.gz"):
            sid_data.export("chimeric/osteosarc-ont.sam", tmp_path / name)
    with pytest.raises(FileNotFoundError):
        sid_data.export("chimeric/osteosarc-ont.sam", tmp_path / "missing" / "reads.sam")


def test_an_unreachable_bundle_fails_once_with_a_hint(tmp_path, monkeypatch):
    import osteosarc
    calls = []

    def unreachable(*args, **kwargs):
        calls.append(args)
        raise OSError("network is unreachable")

    monkeypatch.setattr(osteosarc, "fetch_bundle", unreachable)
    cache = str(tmp_path / "empty-cache")
    for _ in range(2):
        with pytest.raises(sid_data.SidDataUnavailable, match="OSTEOSARC_CACHE"):
            sid_data.bundle(cache)
    assert len(calls) == 1
    sid_data._failures.pop(cache)


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
    dataset = SimpleNamespace(id="snapshot", manifest={"sources": {}}, file=lambda u: asset,
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
