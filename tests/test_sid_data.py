"""Isovar's packaged Sid test reads: openvax-v1's members, record for record, which osteosarc makes again."""

import json
from hashlib import sha256
from pathlib import Path
import shutil
from types import SimpleNamespace

import pysam
import pytest

from isovar import sid_data

DATA = Path(__file__).parent / "data"

# Canonical JSON fingerprints captured before removing any embedded SAM text.
# These pin every decoded value, including order, repeats, headers and metadata.
FIXTURE_DIGESTS = {
    "fusions/coding-corpus/ATP5MG--KMT2A.input.json.gz": "6c4974ae83898705f647f15f2dcdc9e86a913b9baf5577d2326720620df938a1",
    "fusions/corpus/FOXO3--STRADA-CCDC47-T1.input.json.gz": "3136a2d73d9c948f75c388d561b3652239e27e37a1fd80f92d544e8759c6c184",
    "fusions/corpus/FOXO3--STRADA-CCDC47-T2.input.json.gz": "b606792c6d3decd19340ca3fe0aee25fd7c39f8a8743492fec8ecbc714a527bd",
    "fusions/corpus/PARD3B--CDKN2B-AS1-CDKN2B-T1.input.json.gz": "e3a430065c59fadba05b65ffb8b1ff0cb50e35a5dffc2e4d1e12f47dd3927c85",
    "fusions/corpus/TPST1--CRCP-T1.input.json.gz": "e7bde917aeb47d2f8ff76e12a6348dcf4a4347f0357172c4cb417bc95fe15353",
    "fusions/corpus/TPST1--CRCP-T2.input.json.gz": "0094505e4591e88ef7e2103219b006e6961db74856373029db0ec0bbdc847efc",
    "fusions/three-fusions/GABBR1--SLC29A1.ONT-T1.json.gz": "b1f2a3561e22b9159b2fc0b282d70e75c7045fd56a15ad953a1f1cead74566b8",
    "fusions/three-fusions/GABBR1--SLC29A1.ONT-T2.json.gz": "c59fb12794d50d119127570d8b9ab35bfaf5b3953f477e6915188874804f0131",
    "fusions/three-fusions/GABBR1--SLC29A1.ONT-T3.json.gz": "2a65812d89a3f5e47542c891c7d0647462ef6dbc9efe8f845b14eb7b6499bc0b",
    "fusions/three-fusions/GABBR1--SLC29A1.PacBio-T1.json.gz": "6641333dd6e6b248f7bc87323195fb568d9ca8565ee5ca3e7c14bfe3a37f1428",
    "fusions/three-fusions/OTUD7A--FMN1.ONT-T1.json.gz": "3964a67fd55cd6aa068c803acc9362ce6920969d16aee714ac7ecad0f4bc408a",
    "fusions/three-fusions/OTUD7A--FMN1.ONT-T2.json.gz": "0963a0144580bc125af81526b9869776d15b540258478fc6ee2f9e3f97d3bd21",
    "fusions/three-fusions/OTUD7A--FMN1.ONT-T3.json.gz": "3964a67fd55cd6aa068c803acc9362ce6920969d16aee714ac7ecad0f4bc408a",
    "fusions/three-fusions/OTUD7A--FMN1.PacBio-T1.json.gz": "3964a67fd55cd6aa068c803acc9362ce6920969d16aee714ac7ecad0f4bc408a",
    "fusions/three-fusions/PARD3B--CDKN2B-AS1-CDKN2B.ONT-T1.json.gz": "9786df7ab40ae2c429f767b5c2d77f7bf773d866b51fb04d98a069abcb08b40a",
    "fusions/three-fusions/PARD3B--CDKN2B-AS1-CDKN2B.ONT-T2.json.gz": "032ec0b3dcddf764316e97a9f7ca8304b5abc23d00912e9d4f52ee71860c55ee",
    "fusions/three-fusions/PARD3B--CDKN2B-AS1-CDKN2B.ONT-T3.json.gz": "8e4d37a8e0aaeb6cfefdec2e8d144fd080c703b325ef2ffb11deb74ea2133cb2",
    "fusions/three-fusions/PARD3B--CDKN2B-AS1-CDKN2B.PacBio-T1.json.gz": "da5cc08534121475644dc9ab03fe9e913971726342519d276896bee7cf9f9870",
    "osteosarc/figure_comparisons/corpus/context-evidence.json.gz": "0b8bc0d42c5d97d38087f8a99ac3aae1881ee8fa03d6f6c17008bf5c26186a05",
    "osteosarc/figure_comparisons/corpus/dlg5.json.gz": "a3cdbb1dfa278ec1548e8ba8bfee3a7c6c3f89d043aa05d3d56240db31c80b7d",
    "osteosarc/figure_comparisons/corpus/extended-footprints.json.gz": "da721cb27c228f0517806b652cd3c9a79198419ab4619aa2a9d911b968f9861e",
    "osteosarc/figure_comparisons/corpus/insertion-boundaries.json.gz": "0556db50a1414f988c4dd11d81bd4445b2a9305c4734fd800dc76b696067f682",
    "osteosarc/figure_comparisons/corpus/nr2f2-evidence.json.gz": "63bb29726d4a2b7ec75c5597cba97b9b84e0cf88f52968c41ae89b8e38328f0d",
    "osteosarc/figure_comparisons/corpus/nr2f2-libraries.json.gz": "ced1a4ebcd1121d0db5feb1901b3245779d3945ed6f40b7cd0674e9f570f3906",
    "osteosarc/figure_comparisons/corpus/rna-footprints.json.gz": "d69d61c436e7a6790cfc0e9397e07d9c4f063e861e2692a795376e5cad3b462a"
}


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
    fixtures, decoded = {}, {}
    for name in selections:
        path, _, selector = name.partition("#")
        if selector.startswith("/"):
            if path not in decoded:
                data = sid_data.read_fixture(DATA / path)
                fingerprint = sha256(json.dumps(data, sort_keys=True, separators=(",", ":")).encode()).hexdigest()
                assert fingerprint == FIXTURE_DIGESTS[path]
                decoded[path] = tmp_path / ("%d.json" % len(decoded))
                decoded[path].write_text(json.dumps(data))
            fixtures[sid_data.PREFIX + name] = {"json": str(decoded[path]), "pointer": selector}
        else:
            # The read-end regression consumes each selected record directly.
            single = tmp_path / ("%d.sam" % len(fixtures))
            sid_data.export(name, single)
            fixtures[sid_data.PREFIX + name] = str(single)
    assert set(decoded) == set(FIXTURE_DIGESTS)
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


def test_a_missing_or_damaged_packaged_bundle_fails_with_a_hint(tmp_path, monkeypatch):
    packaged = sid_data.PACKAGED
    monkeypatch.setattr(sid_data, "PACKAGED", tmp_path / "missing")
    with pytest.raises(sid_data.SidDataUnavailable, match="python -m isovar.sid_data build"):
        sid_data.bundle()
    changed = tmp_path / "changed"
    changed.mkdir()
    (changed / "manifest.json").write_text("{}")
    monkeypatch.setattr(sid_data, "PACKAGED", changed)
    with pytest.raises(sid_data.SidDataUnavailable, match="reinstall Isovar"):
        sid_data.bundle()
    # A file the manifest doesn't list, as macOS Finder leaves in a checkout.
    damaged = tmp_path / "damaged"
    shutil.copytree(packaged, damaged)
    (damaged / ".DS_Store").write_bytes(b"")
    monkeypatch.setattr(sid_data, "PACKAGED", damaged)
    with pytest.raises(sid_data.SidDataUnavailable, match="unlisted"):
        sid_data.members()


def test_osteosarc_makes_the_packaged_bundle_again(capsys):
    # From the packaged bundle's own recipe and records, offline. CI also makes
    # it again from openvax-v1 (python -m isovar.sid_data check).
    assert sid_data.main(["check", "--from", str(sid_data.PACKAGED)]) is None
    assert "makes Isovar's packaged bundle again" in capsys.readouterr().out


def test_acquiring_again_needs_the_recipes_snapshot(tmp_path):
    with pytest.raises(sid_data.SidDataUnavailable, match="snapshot 2026-09-25"):
        sid_data.build(tmp_path / "bundle", sid_data.PACKAGED, sid=True, cache=str(tmp_path / "empty-cache"))
    assert not (tmp_path / "bundle").exists()


def test_isovars_recipe_keeps_its_members_and_what_they_use():
    isovar_target, other_target = "fixture:isovar/reads.sam", "KRAS-chr12-25245350"
    shared = dict(
        schema_version=1, id="openvax-v2", kind="shared", snapshot={"id": "snapshot"},
        aliases={"old-reads": isovar_target, "KRAS-old": other_target},
        targets={isovar_target: {"kind": "fixture"}, other_target: {"kind": "small_variant"}},
        sources={"rna": {"sample": "T1"}, "wgs": {"sample": "T2"}},
        members={"isovar/reads.sam": dict(source="rna", target=isovar_target),
                 "varcode/reads.sam": dict(source="wgs", target=other_target),
                 "topiary/isovar/reads.sam": dict(source="wgs", target=other_target)})
    mine = sid_data.recipe(shared)
    assert mine == dict(shared, id="isovar", aliases={"old-reads": isovar_target},
                        targets={isovar_target: {"kind": "fixture"}}, sources={"rna": {"sample": "T1"}},
                        members={"isovar/reads.sam": dict(source="rna", target=isovar_target)})
    assert sid_data.recipe(mine) == mine
    assert "aliases" not in sid_data.recipe({k: v for k, v in shared.items() if k != "aliases"})


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


def test_record_references_preserve_order_repeats_and_metadata(monkeypatch, tmp_path):
    import socket

    def no_network(*args, **kwargs):
        raise AssertionError("Fixture lookup must stay offline")

    monkeypatch.setattr(socket, "create_connection", no_network)
    name = "fusions/coding-corpus/ATP5MG--KMT2A.input.json.gz#/original_records"
    lookup = sid_data.records(name, cache=str(tmp_path))
    first, second = lookup
    reference = lambda digest: {"sid_member": name, "sam_sha256": digest}
    encoded = [{"sam": reference(second), "partner_sam": reference(first),
                "minimum_base_quality": 27, "source_reverse_complement": True},
               reference(second), {"sam_sha256": first}, None]
    saved = json.dumps(encoded)
    assert sid_data.restore_records(encoded, cache=str(tmp_path)) == [
        {"sam": lookup[second], "partner_sam": lookup[first],
         "minimum_base_quality": 27, "source_reverse_complement": True},
        lookup[second], {"sam_sha256": first}, None]
    assert json.dumps(encoded) == saved


def test_record_references_reject_missing_and_wrong_member_digests():
    name = "fusions/coding-corpus/ATP5MG--KMT2A.input.json.gz#/original_records"
    digest = next(iter(sid_data.records(name)))
    for member, checksum in [("no-such-member", digest), (name, "0" * 64),
                             ("chimeric/osteosarc-ont.sam", digest)]:
        with pytest.raises(KeyError):
            sid_data.restore_records({"sid_member": member, "sam_sha256": checksum})


@pytest.mark.parametrize("reference", [
    {"sid_member": "member"}, {"sid_member": [], "sam_sha256": "digest"},
    {"sid_member": "member", "sam_sha256": []},
    {"sid_member": "member", "sam_sha256": "digest", "extra": True},
])
def test_malformed_record_reference_is_not_silently_accepted(reference):
    with pytest.raises(ValueError, match="Invalid Sid record reference"):
        sid_data.restore_records(reference)


def test_record_lookup_rejects_changed_cached_bam(tmp_path, monkeypatch):
    name = "chimeric/osteosarc-ont.sam"
    source = sid_data._file(name, "bam", str(tmp_path))
    changed = tmp_path / "changed.bam"
    with pysam.AlignmentFile(source) as reader, pysam.AlignmentFile(changed, "wb", template=reader) as writer:
        for read in reader:
            read.mapping_quality = (read.mapping_quality + 1) % 255
            writer.write(read)
    monkeypatch.setattr(sid_data, "_file", lambda *args: changed)
    with pytest.raises(sid_data.SidDataUnavailable, match="Exported Sid records differ"):
        sid_data.records(name)


def test_converted_fixtures_do_not_embed_sam_records():
    def check(value):
        if isinstance(value, str):
            assert value.count("\t") < 10
        elif isinstance(value, list):
            for child in value:
                check(child)
        elif isinstance(value, dict):
            for child in value.values():
                check(child)

    for filename in FIXTURE_DIGESTS:
        check(sid_data.read_json(DATA / filename))
    assert not (DATA / "read_ends/osteosarc.sam").exists()
