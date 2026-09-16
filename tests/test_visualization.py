"""Figure geometry and evidence stay faithful to the existing protein pipeline."""

from copy import deepcopy
import json
from pathlib import Path
import re
import subprocess
import sys
from types import SimpleNamespace

import pytest
from varcode import Variant

from isovar import ProteinSequenceCreator
from isovar.default_parameters import PLOT_DPI, PLOT_MAX_ROWS, PLOT_VIEW
from isovar.read_evidence import ReadEvidence
from isovar.visualization import (
    _creator_settings, _genomic_projection, _witness_data, collect_visualization_data,
    plot_variant_evidence, save_variant_figures, timestamped_run_directory, variant_directory_name,
)
from .test_translation_regressions import _context, _read


@pytest.fixture(params=["+", "-"])
def evidence_data(request, monkeypatch):
    strand = request.param
    variant = Variant("1", 1000, "G", "C", genome="GRCh38")
    prefix, suffix = "AAA" * 30, "CCC" * 30
    context = _context(variant, strand, prefix, suffix)
    # Unknown path metadata remains unclassified until the normal pipeline.
    context.transcripts = (SimpleNamespace(id="TX1", gene_name="EXAMPLE", strand=strand,
                                          exon_intervals=((1, 2000),)),)
    monkeypatch.setattr("isovar.protein_sequence_creator.reference_contexts_for_variant", lambda *a, **kw: [context])
    alt = "C" if strand == "+" else "G"
    reads = [_read(prefix, alt, suffix[:30], "left%d" % i, strand) for i in range(3)]
    reads += [_read(prefix[-30:], alt, suffix, "right%d" % i, strand) for i in range(3)]
    evidence = ReadEvidence(1000, "G", "C", [], reads, [])
    data = collect_visualization_data(variant, evidence, compare_assembly=True)
    return data, variant, evidence


def test_comparison_is_exactly_the_existing_pipeline(evidence_data):
    data, variant, evidence = evidence_data
    for mode in data["modes"]:
        p = ProteinSequenceCreator(variant_sequence_assembly=mode["assembly"]).sorted_protein_sequences_for_variant(
            variant, evidence)[0]
        assert mode["protein"]["amino_acids"] == p.amino_acids
        assert mode["protein"]["templates"] == p.num_supporting_fragments
    a, b = [m["settings"] for m in data["modes"]]
    assert {k for k in a if a[k] != b[k]} == {"variant_sequence_assembly"}
    assert len(data["modes"][0]["protein"]["amino_acids"]) > len(data["modes"][1]["protein"]["amino_acids"])
    json.dumps(data)  # No NumPy scalars or opaque sequence/transcript objects.


def test_oriented_span_coverage_is_exact(evidence_data):
    data, _, _ = evidence_data
    for mode in data["modes"]:
        w = mode["protein"]["witness"]
        expected = [sum(s["observations"] for s in w["spans"] if s["start"] <= x < s["end"])
                    for x in range(w["start"], w["end"])]
        assert w["coverage"] == expected
        assert len(w["cdna"]) == w["end"] - w["start"]
        assert sum(s["observations"] for s in w["spans"]) == w["observations"]
        assert "left0" not in json.dumps(data)
    assert data["modes"][0]["protein"]["witness"]["spanning_observations"] == 0


@pytest.mark.parametrize("view", ["all", "protein", "assembly", "transcripts"])
def test_all_views_render_opaque_and_do_not_mutate_evidence(evidence_data, view):
    pytest.importorskip("matplotlib")
    from matplotlib.backends.backend_agg import FigureCanvasAgg

    data, _, _ = evidence_data
    before = deepcopy(data)
    figure = plot_variant_evidence(data, view=view, max_rows=2)
    FigureCanvasAgg(figure).draw()
    assert figure.get_facecolor() == (1, 1, 1, 1)
    assert data == before
    assert all(ax.get_facecolor() == (1, 1, 1, 1) for ax in figure.axes)
    if view in {"all", "assembly"}:
        offset = 1 if view == "all" else 0
        assert figure.axes[offset].get_xlim() == figure.axes[offset + 1].get_xlim()


def test_no_protein_is_explicit_not_fabricated(monkeypatch):
    pytest.importorskip("matplotlib")
    from matplotlib.backends.backend_agg import FigureCanvasAgg

    variant = Variant("1", 1000, "G", "C", genome="GRCh38")
    evidence = ReadEvidence(1000, "G", "C", [], [], [])
    data = collect_visualization_data(variant, evidence, compare_assembly=True)
    assert all(m["protein"] is None for m in data["modes"])
    figure = plot_variant_evidence(data)
    FigureCanvasAgg(figure).draw()
    assert "No translated protein" in [t.get_text() for ax in figure.axes for t in ax.texts]


def test_compressed_axis_preserves_exons_and_is_monotonic():
    data = dict(variant=dict(interval=[110, 111]),
                transcripts=[dict(exons=[(100, 150), (1000, 1050)])],
                modes=[dict(protein=dict(witness=dict(genomic_blocks=[(100, 150), (1000, 1050)])))])
    lo, hi, segments, project = _genomic_projection(data)
    assert project(150) - project(100) == 50
    assert project(1000) - project(150) == 25
    assert project(1050) - project(1000) == 50
    values = [project(x) for x in range(lo, hi + 1)]
    assert values == sorted(values)
    with pytest.raises(ValueError):
        project(hi + 1)


def test_svg_png_json_and_timestamp_layout(evidence_data, tmp_path):
    pytest.importorskip("matplotlib")
    from PIL import Image

    data, _, _ = evidence_data
    run = timestamped_run_directory(tmp_path)
    assert re.fullmatch(r"\d{4}-\d{2}-\d{2}_\d{2}-\d{2}-\d{2}-\d{6}Z", run.name)
    assert timestamped_run_directory(tmp_path) != run
    directory = save_variant_figures(data, run, view="protein", dpi=72)
    assert directory.parent == run
    svg = (directory / "protein.svg").read_text()
    assert "<text" in svg and "Assembly on" in svg and "#ffffff" in svg
    with Image.open(directory / "protein.png") as image:
        assert image.convert("RGBA").getpixel((0, 0)) == (255, 255, 255, 255)
    saved = json.loads((directory / "evidence.json").read_text())
    assert saved["modes"] == data["modes"]
    with pytest.raises(FileExistsError):
        save_variant_figures(data, run, view="protein", dpi=72)


def test_plotting_is_a_lazy_optional_dependency():
    result = subprocess.run([sys.executable, "-c",
                             "import isovar.visualization, sys; assert 'matplotlib' not in sys.modules"],
                            capture_output=True, text=True)
    assert result.returncode == 0, result.stderr


def test_cli_defaults_match_plotting_and_protein_apis():
    from isovar.cli.isovar_plot import parser
    from isovar.cli.translation_args import protein_sequence_creator_kwargs_from_args

    args = parser.parse_args(["--bam", "input.bam", "--variant", "1", "1000", "G", "C", "--genome", "GRCh38"])
    assert (args.view, args.dpi, args.max_rows) == (PLOT_VIEW, PLOT_DPI, PLOT_MAX_ROWS)
    assert _creator_settings(ProteinSequenceCreator(**protein_sequence_creator_kwargs_from_args(args))) == \
        _creator_settings(ProteinSequenceCreator())


@pytest.mark.parametrize("options", [{"max_rows": 1}, {"max_rows": True}, {"view": "invalid"}])
def test_invalid_render_options(evidence_data, options):
    with pytest.raises(ValueError):
        plot_variant_evidence(evidence_data[0], **options)


def test_cli_rejects_multiple_variants_before_opening_bam(monkeypatch):
    from isovar.cli import isovar_plot

    monkeypatch.setattr(isovar_plot, "variant_collection_from_args", lambda args: [1, 2])
    with pytest.raises(SystemExit) as error:
        isovar_plot.run(["--bam", "does-not-exist.bam"])
    assert error.value.code == 2


def test_long_alleles_have_bounded_distinct_directory_names(evidence_data):
    data = deepcopy(evidence_data[0])
    data["variant"]["alt"] = "A" * 1000
    a = variant_directory_name(data)
    data["variant"]["alt"] += "T"
    b = variant_directory_name(data)
    assert a != b and max(len(a), len(b)) <= 120


def test_junctions_outside_the_witness_are_not_displayed():
    from isovar.allele_read import AlleleRead
    from isovar.variant_sequence import VariantSequence

    read = AlleleRead("A" * 10, "C", "T" * 10, "r", reference_blocks=(
        (0, 6, 90, 96), (6, 10, 200, 204), (10, 21, 300, 311)),
        splice_junctions=((96, 200), (204, 300)))
    sequence = VariantSequence("AAA", "C", "TTT", [read])
    orf = SimpleNamespace(cdna_sequence="AAACTTT", offset_to_first_complete_codon=0,
                          variant_cdna_interval_start=3, variant_cdna_interval_end=4)
    translation = SimpleNamespace(untrimmed_variant_sequence=sequence, variant_orf=orf,
                                  reference_context=SimpleNamespace(strand="+", transcripts=[]))
    w = _witness_data(translation)
    assert w["genomic_blocks"] == [(201, 204), (300, 304)]
    assert w["junctions"] == [dict(start=204, end=300, observations=1)]


def test_deletion_boundary_and_long_protein_render(evidence_data):
    pytest.importorskip("matplotlib")
    from matplotlib.backends.backend_agg import FigureCanvasAgg

    data = deepcopy(evidence_data[0])
    for mode in data["modes"]:
        mode["protein"]["amino_acids"] = "A" * 150
        mode["protein"]["mutation_start"] = mode["protein"]["mutation_end"] = 20
    figure = plot_variant_evidence(data, view="protein")
    FigureCanvasAgg(figure).draw()
    assert any("sequence in evidence.json" in t.get_text() for t in figure.axes[0].texts)


def test_cli_writes_a_timestamped_report(evidence_data, monkeypatch, tmp_path):
    pytest.importorskip("matplotlib")
    from contextlib import nullcontext
    from isovar.cli import isovar_plot

    data, variant, evidence = evidence_data
    monkeypatch.setattr(isovar_plot, "variant_collection_from_args", lambda args: [variant])
    monkeypatch.setattr(isovar_plot, "alignment_file_from_args", lambda args: nullcontext(None))
    monkeypatch.setattr(isovar_plot, "read_collector_from_args", lambda args: SimpleNamespace(
        read_evidence_for_variant=lambda *a: evidence, min_mapping_quality=1,
        use_duplicate_reads=False, use_secondary_alignments=True, use_soft_clipped_bases=False,
        merge_overlapping_fragments=True))
    isovar_plot.run(["--bam", "rna.bam", "--compare-assembly", "--view", "protein",
                     "--dpi", "72", "--output-dir", str(tmp_path)])
    reports = list(tmp_path.glob("*/*/evidence.json"))
    assert len(reports) == 1
    assert json.loads(reports[0].read_text())["modes"] == data["modes"]
