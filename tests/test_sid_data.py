"""Isovar's Sid test reads come from openvax-v1, record for record."""

from collections import Counter
import gzip
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


def _embedded_lines(value):
    """SAM lines inside a JSON fixture value: strings, lists, {digest: line} maps and sam/partner_sam fields."""
    if isinstance(value, str):
        if len(value.split("\t")) >= 11:
            yield value
    elif isinstance(value, list):
        for item in value:
            yield from _embedded_lines(item)
    elif isinstance(value, dict):
        for key, item in value.items():
            if key in ("sam", "partner_sam") or isinstance(item, (list, dict)) or len(str(key)) == 64:
                yield from _embedded_lines(item)


def _pointer(document, pointer):
    for part in pointer.lstrip("/").split("/"):
        if part:
            part = part.replace("~1", "/").replace("~0", "~")
            document = document[int(part)] if isinstance(document, list) else document[part]
    return document


def test_every_read_file_is_exported_under_its_original_name(members):
    files = [name for name in members if "#" not in name]
    assert len(files) == 124
    for name in files:
        path = sid_data.path(name)
        assert path.name == Path(name).name
        if name.endswith(".bam"):
            assert Path(str(path) + ".bai").exists()
        _records(path)
    manifest = json.loads((DATA / "osteosarc/manifest.json").read_text())
    for data in manifest["datasets"].values():
        _, records = _records(sid_data.path("osteosarc/" + data["file"]))
        assert len(records) == data["fixture_records"]


def test_embedded_and_selected_reads_match_the_bundle(members, tmp_path):
    from osteosarc import export_bundle
    selections = [name for name in members if "#" in name]
    assert len(selections) == 221
    exported = export_bundle(sid_data.bundle(), tmp_path, format="bam",
                             members=[sid_data.PREFIX + name for name in selections])
    for name in selections:
        header, expected = _records(exported[sid_data.PREFIX + name])
        path, _, selector = name.partition("#")
        if selector.startswith("/"):
            opener = gzip.open if path.endswith(".gz") else open
            with opener(DATA / path, "rt") as handle:
                lines = list(_embedded_lines(_pointer(json.load(handle), selector)))
        else:
            # One read of a SAM file; ":" in its name is "-" in the member name.
            _, records = _records(DATA / path)
            lines = [r for r in records if r.split("\t")[0] in (selector, selector.replace("-", ":"))]
        # Compare through one header, so text details such as float formatting agree.
        observed = [pysam.AlignedSegment.fromstring(line, header).to_string() for line in lines]
        assert Counter(observed) == Counter(expected), name


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


def test_damaged_exports_are_detected(tmp_path):
    for name, text in (("a.sam", "one"), ("b/c.sam", "two")):
        (tmp_path / name).parent.mkdir(parents=True, exist_ok=True)
        (tmp_path / name).write_text(text)
    (tmp_path / "digests.json").write_text(json.dumps(sid_data._digests(tmp_path, ["a.sam", "b/c.sam"])))
    sid_data._verify_exports(tmp_path, {"a.sam", "b/c.sam"})
    (tmp_path / "b/c.sam").write_text("edited")
    with pytest.raises(sid_data.SidDataUnavailable, match="delete it"):
        sid_data._verify_exports(tmp_path, {"a.sam", "b/c.sam"})
    (tmp_path / "b/c.sam").unlink()
    with pytest.raises(sid_data.SidDataUnavailable):
        sid_data._verify_exports(tmp_path, {"a.sam", "b/c.sam"})


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
