"""The selected figure captions must remain independently reproducible."""

import json

from examples import osteosarc_assembly_figures


def test_osteosarc_figure_comparisons_validate_without_changing_reads(tmp_path, monkeypatch):
    # Rendering itself is tested separately; keep this scientific regression
    # usable even when the optional Matplotlib dependency is not installed.
    from isovar.visualization import variant_directory_name

    def record_data(data, output):
        directory = output / variant_directory_name(data)
        directory.mkdir()
        (directory / "evidence.json").write_text(json.dumps(data))
        return directory

    monkeypatch.setattr(osteosarc_assembly_figures, "save_variant_figures", record_data)
    output = osteosarc_assembly_figures.generate(tmp_path)
    manifest = json.loads((output / "manifest.json").read_text())
    assert [(e["on_length"], e["off_length"], e["on_windows"], e["off_windows"])
            for e in manifest["examples"]] == [(49, 40, 25, 16), (49, 44, 25, 20), (49, 37, 25, 13)]
    for path in output.glob("*/evidence.json"):
        data = json.loads(path.read_text())
        assert data["modes"][0]["protein"]["witness"]["spanning_observations"] == 0
        assert data["provenance"]["source_bam_url"].startswith("https://")
        assert data["provenance"]["nonfocal_witness_indels"] == 0
        assert all(c["matches_expected"] for v in data["provenance"]["independent_validation"] for c in v["checks"])
