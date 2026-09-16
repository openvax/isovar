"""The selected figure captions must remain independently reproducible."""

import json
from hashlib import sha256

import pysam

from examples import osteosarc_assembly_figures


def record_data(data, output):
    from isovar.visualization import variant_directory_name

    directory = output / variant_directory_name(data)
    directory.mkdir(parents=True)
    (directory / "evidence.json").write_text(json.dumps(data))
    return directory


def test_osteosarc_figure_comparisons_validate_without_changing_reads(tmp_path, monkeypatch):
    # Rendering itself is tested separately; keep this scientific regression
    # usable even when the optional Matplotlib dependency is not installed.
    monkeypatch.setattr(osteosarc_assembly_figures, "save_variant_figures", record_data)
    output = osteosarc_assembly_figures.generate(tmp_path, case_ids=osteosarc_assembly_figures.CASE_IDS)
    manifest = json.loads((output / "manifest.json").read_text())
    assert [(e["on_length"], e["off_length"], e["on_windows"], e["off_windows"])
            for e in manifest["examples"]] == [(49, 40, 25, 16), (49, 44, 25, 20), (49, 37, 25, 13)]
    for path in output.glob("*/evidence.json"):
        data = json.loads(path.read_text())
        assert all(t["name"] and t["name"] != t["id"] for t in data["transcripts"])
        junctions = {(j["start"], j["end"]) for j in data["modes"][0]["protein"]["witness"]["junctions"]}
        assert junctions
        for transcript in data["transcripts"]:
            if not transcript["rna_supported"]:
                continue
            introns = {(a[1], b[0]) for a, b in zip(transcript["exons"], transcript["exons"][1:])}
            assert junctions <= introns
        assert data["modes"][0]["protein"]["witness"]["spanning_observations"] == 0
        assert data["provenance"]["source_bam_url"].startswith("https://")
        assert data["provenance"]["nonfocal_witness_indels"] == 0
        assert all(c["matches_expected"] for v in data["provenance"]["independent_validation"] for c in v["checks"])


def test_extended_examples_separate_predictions_from_rna_and_keep_ambiguity(tmp_path, monkeypatch):
    monkeypatch.setattr(osteosarc_assembly_figures, "save_variant_figures", record_data)
    output = osteosarc_assembly_figures.generate(tmp_path, case_ids=osteosarc_assembly_figures.ADDITIONAL_CASE_IDS)
    cases = {d["provenance"]["case_id"].split("-")[1]: d
             for path in output.rglob("evidence.json") for d in [json.loads(path.read_text())]}
    map2 = cases["MAP2"]
    assert map2["counts"]["alt"]["observations"] == 1
    assert all(m["protein"] is None for m in map2["modes"])
    assert len([p for p in map2["reference_predictions"] if p["protein"]]) == 2
    assert all(not t["rna_supported"] for t in map2["transcripts"])
    nav2 = cases["NAV2"]
    assert nav2["modes"][0]["protein"]["transcript_ids"] == ["ENST00000527559", "ENST00000540292"]
    assert len([p for p in nav2["reference_predictions"] if p["protein"]]) == 7
    assert sum(t["rna_supported"] for t in nav2["transcripts"]) == 2
    assert len(nav2["modes"][0]["protein"]["amino_acids"]) == 47
    assert nav2["modes"][0]["protein"]["amino_acids"].endswith("W")
    assert nav2["modes"][0]["protein"]["witness"]["junctions"] == []
    assert {p["protein"]["amino_acids"][-3:] for p in nav2["reference_predictions"] if p["protein"]} == {"WLR", "WVN"}
    ntf3 = cases["NTF3"]
    for mode in ntf3["modes"]:
        p = mode["protein"]
        assert p["amino_acids"][p["mutation_start"]] == "S"
        witness = p["witness"]
        offset = -witness["start"]
        assert witness["cdna"][offset:offset + 2] == "GT"  # Original genomic AG>GT compound, not just a matching amino acid.
    for prediction in ntf3["reference_predictions"]:
        p = prediction["protein"]
        assert p["amino_acids"][p["mutation_start"]] == "R"
    assert {p["description"] for p in ntf3["reference_predictions"]} == {"p.K56R", "p.K69R"}
    for data in cases.values():
        assert all(v["status"] in {"ok", "no_protein"} for v in data["provenance"]["independent_validation"])
        assert data["provenance"]["nonfocal_witness_indels"] == 0
    assert all(not c["matches_expected"] for v in ntf3["provenance"]["independent_validation"] for c in v["checks"])


def test_t1_long_short_comparison_uses_complete_original_regions(tmp_path, monkeypatch):
    monkeypatch.setattr(osteosarc_assembly_figures, "save_variant_figures", record_data)
    corpus = osteosarc_assembly_figures.PAIR_CORPUS
    for case in json.loads((corpus / "manifest.json").read_text())["cases"]:
        with pysam.AlignmentFile(corpus / case["primary_bam"]) as bam:
            reads = list(bam)
        assert [sha256(r.to_string().encode()).hexdigest() for r in reads] == case["retained_record_sha256"]
        assert all(not r.flag & (256 | 1024 | 2048) for r in reads)
        assert len(case["original_record_sha256"]) == case["acquisition"]["region_records"]
    output = osteosarc_assembly_figures.generate(tmp_path, case_ids=osteosarc_assembly_figures.PAIR_CASE_IDS)
    datasets = {d["provenance"]["case_id"]: d
                for path in output.rglob("evidence.json") for d in [json.loads(path.read_text())]}
    long = datasets["PIP5K1A-T1-ONT"]
    short = datasets["PIP5K1A-T1-Illumina"]
    assert long["counts"]["alt"]["observations"] == 4
    assert short["counts"]["alt"]["observations"] == 1
    assert all(m["protein"] is None for m in short["modes"])
    assert {m["protein"]["amino_acids"] for m in long["modes"]} == {
        "VFKKIPLKPSPSKKFRSGSSFSRRAAPVATPALLTSHRSLGNTRHK"}
    assert all(c["matches_expected"] for v in long["provenance"]["independent_validation"] for c in v["checks"])
    for a, b in zip(long["modes"], short["modes"]):
        assert a["settings"] == b["settings"]
