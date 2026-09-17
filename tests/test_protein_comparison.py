"""Complete, source-labelled alternatives without changing default selection."""

from copy import deepcopy
import json

import pytest

from isovar import ProteinSequenceCreator
from isovar.protein_comparison import comparison_rows, protein_groups, protein_comparison_figures, save_protein_comparison
from isovar.visualization import collect_visualization_data
from .test_visualization import evidence_data as evidence_data


def uncapped(evidence_data):
    _, variant, evidence = evidence_data
    return collect_visualization_data(variant, evidence, compare_assembly=True,
                                      creator_kwargs=dict(max_protein_sequences_per_variant=None))


def test_all_returned_proteins_and_frame_contexts_match_pipeline(evidence_data):
    data = uncapped(evidence_data)
    _, variant, evidence = evidence_data
    for mode in data["modes"]:
        expected = ProteinSequenceCreator(variant_sequence_assembly=mode["assembly"],
                                          max_protein_sequences_per_variant=None).sorted_protein_sequences_for_variant(variant, evidence)
        assert [p["amino_acids"] for p in mode["proteins"]] == [p.amino_acids for p in expected]
        assert mode["protein"] == next(iter(mode["proteins"]), None)
        for protein, original in zip(mode["proteins"], expected):
            assert len(protein["frames"]) == len(original.translations)
            for frame, translation in zip(protein["frames"], original.translations):
                orf = translation.variant_orf
                assert frame["cdna"] == orf.cdna_sequence
                assert frame["variant_codon_phase"] == (orf.variant_cdna_interval_start - orf.offset_to_first_complete_codon) % 3
    assert all(len(m["proteins"]) <= 1 for m in evidence_data[0]["modes"])
    selected = {tid for m in data["modes"] if m["protein"] for tid in m["protein"]["transcript_ids"]}
    assert {t["id"] for t in data["transcripts"] if t["rna_supported"]} == selected
    json.dumps(data)


def test_sources_and_assembly_modes_never_pool_template_counts(evidence_data):
    data = uncapped(evidence_data)
    products = [dict(source=s, label=s + " / ONT", visualization=deepcopy(data)) for s in ("T1", "T2")]
    rows = comparison_rows(products)
    for product in products:
        rna = [r for r in rows if r["source"] == product["source"]]
        assert [r["assembly"] for r in rna] == sorted(r["assembly"] for r in rna)
        for mode in data["modes"]:
            assert {rank for row in rna if row["assembly"] == mode["assembly"] for rank in row["covered_ranks"]} == set(range(1, len(mode["proteins"]) + 1))
        for row in rna:
            mode = next(m for m in data["modes"] if m["assembly"] == row["assembly"])
            assert row["protein"] == mode["proteins"][row["rank"] - 1]


def test_partial_context_keeps_multiple_compatible_alternatives_and_frames():
    def p(sequence, phase=0, stop=False):
        return dict(amino_acids=sequence, mutation_start=0, frameshift=False,
                    ends_with_stop_codon=stop, frames=[dict(strand="+", variant_codon_phase=phase, transcript_ids=["TX1"])])
    proteins = [p("AC"), p("ACD"), p("ACE"), p("AC", phase=1), p("AC", stop=True)]
    assert protein_groups(proteins) == [(2, [1, 2]), (3, [1, 3]), (4, [4]), (5, [1, 5])]
    assert protein_groups([p("AC"), p("CD")]) == [(1, [1]), (2, [2])]


def test_short_comparison_pages_keep_notes_clear_of_the_axis(evidence_data):
    from matplotlib.backends.backend_agg import FigureCanvasAgg
    products = [dict(source="T1", label="T1 / ONT", visualization=uncapped(evidence_data))]
    for _, figure in protein_comparison_figures(products, rows_per_page=3):
        canvas = FigureCanvasAgg(figure)
        canvas.draw()
        renderer = canvas.get_renderer()
        xlabel = figure.axes[0].xaxis.label.get_window_extent(renderer)
        for note in figure.texts[1:]:
            assert not xlabel.overlaps(note.get_window_extent(renderer))
        figure.clear()


def test_capped_or_incompatible_sources_cannot_claim_all_alternatives(evidence_data):
    data = uncapped(evidence_data)
    capped = dict(source="T1", label="T1", visualization=evidence_data[0])
    with pytest.raises(ValueError, match="uncapped"):
        comparison_rows([capped])
    a = dict(source="T1", label="T1", visualization=data)
    b = deepcopy(a)
    with pytest.raises(ValueError, match="uniquely"):
        comparison_rows([a, b])
    b["source"] = "T2"
    b["visualization"]["variant"]["start"] += 1
    with pytest.raises(ValueError, match="different variants"):
        comparison_rows([a, b])
    b = deepcopy(a)
    b["source"] = "T2"
    b["visualization"]["transcripts"][0]["exons"] = [[1, 999999]]
    with pytest.raises(ValueError, match="Transcript models differ"):
        comparison_rows([a, b])


def test_cli_all_proteins_is_opt_in_and_keeps_the_top_result(evidence_data, monkeypatch, tmp_path):
    from contextlib import nullcontext
    from isovar.cli import isovar_plot
    from isovar import ReadCollector
    from isovar.default_parameters import PLOT_ALL_PROTEINS
    import isovar.protein_comparison as comparison

    _, variant, evidence = evidence_data
    assert isovar_plot.parser.parse_args(["--bam", "input.bam"]).all_proteins == PLOT_ALL_PROTEINS is False
    monkeypatch.setattr(isovar_plot, "variant_collection_from_args", lambda args: [variant])
    monkeypatch.setattr(isovar_plot, "alignment_file_from_args", lambda args: nullcontext(None))
    monkeypatch.setattr(ReadCollector, "read_evidence_for_variant", lambda *args: evidence)
    captured = {}
    def save(data, *args, **kwargs):
        captured["data"] = data
        return tmp_path
    monkeypatch.setattr(isovar_plot, "save_variant_figures", save)
    monkeypatch.setattr(comparison, "save_protein_comparison", lambda products, *args, **kw: captured.update(products=products))
    isovar_plot.run(["--bam", "input.bam", "--output-dir", str(tmp_path), "--compare-assembly", "--all-proteins", "--sample-label", "T1 / ONT"])
    assert captured["products"][0]["label"] == "T1 / ONT"
    for mode in captured["data"]["modes"]:
        assert mode["settings"]["max_protein_sequences_per_variant"] is None
        assert mode["protein"] == mode["proteins"][0]


def test_no_coverage_no_alternate_and_no_frame_remain_distinct(evidence_data):
    data = uncapped(evidence_data)
    for mode in data["modes"]:
        mode.update(protein=None, proteins=[])
    products = [dict(source="T1", label="T1 / ONT", visualization=data)]
    assert {r["status"] for r in comparison_rows(products) if r["source"] == "T1"} == {"no_translated_protein"}
    data["counts"]["alt"]["templates"] = 0
    data["counts"]["ref"]["templates"] = 1
    assert comparison_rows(products)[0]["status"] == "no_alternate_support"
    for count in data["counts"].values():
        count["templates"] = 0
    assert comparison_rows(products)[0]["status"] == "no_callable_reads"


def test_missing_varcode_prediction_and_pagination_are_not_silent(evidence_data, tmp_path):
    pytest.importorskip("matplotlib")
    from matplotlib.backends.backend_agg import FigureCanvasAgg
    data = uncapped(evidence_data)
    data["reference_predictions"] = [dict(transcript_id="TX1", description="splice region",
                                           protein=None, unavailable_reason="Splicing consequence unresolved")]
    products = [dict(source="T1", label="T1 / ONT", visualization=data)]
    before = deepcopy(products)
    rows = comparison_rows(products)
    assert rows[-1]["status"] == "no_concrete_prediction"
    pages = list(protein_comparison_figures(products, rows_per_page=2))
    assert len(pages) == len([r for r in rows if r["source"] != "Varcode"])
    for _, figure in pages:
        FigureCanvasAgg(figure).draw()
        assert figure.get_facecolor() == (1, 1, 1, 1)
    assert any("No concrete prediction" in t.get_text() for t in pages[-1][1].axes[0].texts)
    directory = save_protein_comparison(products, tmp_path / "comparison", dpi=50, rows_per_page=2)
    assert len(list(directory.glob("*.svg"))) == len(pages)
    assert len(json.loads((directory / "evidence.json").read_text())["rows"]) == len(rows)
    assert products == before


@pytest.mark.parametrize("rows_per_page", [True, 1, 2.5])
def test_invalid_page_size(evidence_data, rows_per_page):
    products = [dict(source="T1", label="T1", visualization=uncapped(evidence_data))]
    with pytest.raises(ValueError):
        list(protein_comparison_figures(products, rows_per_page))
