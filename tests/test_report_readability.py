"""Report labels retain producer identities, outcomes and visible alternatives."""

from copy import deepcopy

import pytest
from varcode import Variant

from isovar.read_evidence import ReadEvidence
from isovar.visualization import collect_visualization_data, plot_variant_evidence
from tests.test_translation_regressions import _read
from tests.test_visualization import evidence_data as evidence_data


def test_focal_counts_union_original_fragments_and_explain_empty_protein():
    variant = Variant('1', 1000, 'G', 'C', genome='GRCh38')
    # One template has discordant allele observations: it belongs to both
    # numerator sets but is only one original template in the denominator.
    ref = [_read('AAA', 'G', 'CCC', 'shared', '+'), _read('AAA', 'G', 'CCC', 'ref-only', '+')]
    alt = [_read('AAA', 'C', 'CCC', 'shared', '+')]
    data = collect_visualization_data(variant, ReadEvidence(1000, 'G', 'C', ref, alt, []),
        creator_kwargs={'min_variant_sequence_coverage': 10}, include_reference_predictions=False)
    assert data['allele_support']['total']['fragments'] == 2
    assert data['allele_support']['alt']['fragments'] == 1
    assert data['allele_support']['ref']['fragments'] == 2
    assert data['modes'][0]['assembly_state'] == 'no_retained_rna'
    assert 'coverage floor 10' in data['modes'][0]['no_protein_reason']
    figure = plot_variant_evidence(data)
    text = '\n'.join(t.get_text() for t in figure.texts + [t for ax in figure.axes for t in ax.texts])
    assert 'ALT 1/2 fragments (REF 2; other 0)' in text
    assert 'coverage floor 10' in text
    assert 'alternate read objects' not in text


def test_transcript_cap_keeps_supported_then_focal_exon_models(evidence_data):
    data = deepcopy(evidence_data[0])
    supported = dict(data['transcripts'][0], id='supported', rna_supported=True)
    overlapping = dict(supported, id='focal', rna_supported=False)
    irrelevant = [dict(supported, id=f'blank-{i}', rna_supported=False, exons=[(1, 10)]) for i in range(25)]
    data['transcripts'] = irrelevant + [overlapping, supported]
    figure = plot_variant_evidence(data, view='transcripts', max_rows=2)
    text = '\n'.join(t.get_text() for t in figure.axes[0].texts)
    assert 'supported' in text and 'focal' in text
    assert 'blank-' not in text
    assert data['transcripts'][0]['id'] == 'blank-0'


@pytest.mark.parametrize('view', ['all', 'protein', 'transcripts'])
def test_allele_counts_do_not_overlap_the_first_panel_title(evidence_data, view):
    from matplotlib.backends.backend_agg import FigureCanvasAgg
    from matplotlib.text import Text

    figure = plot_variant_evidence(evidence_data[0], view=view)
    canvas = FigureCanvasAgg(figure)
    canvas.draw()
    note, = [text for text in figure.texts if text.get_text().startswith('Focal RNA:')]
    axis = figure.axes[0]
    title, = [text for text in axis.get_children()
              if isinstance(text, Text) and text.get_text() == axis.get_title(loc='left')]
    assert title.get_window_extent(canvas.get_renderer()).y1 < note.get_window_extent(canvas.get_renderer()).y0


def test_plot_selects_one_loaded_variant_before_collecting_rna(evidence_data, monkeypatch, tmp_path):
    from contextlib import nullcontext
    from isovar.cli import isovar_plot
    from isovar import ReadCollector

    _, variant, evidence = evidence_data
    monkeypatch.setattr(isovar_plot, 'variant_collection_from_args', lambda args: [object(), variant])
    monkeypatch.setattr(isovar_plot, 'alignment_file_from_args', lambda args: nullcontext(None))
    selected = []
    def collect(self, observed_variant, alignment):
        selected.append(observed_variant)
        return evidence
    monkeypatch.setattr(ReadCollector, 'read_evidence_for_variant', collect)
    exported = []
    def save(data, *args, **kwargs):
        exported.append(data)
        return tmp_path
    monkeypatch.setattr(isovar_plot, 'save_variant_figures', save)
    isovar_plot.run(['--bam', 'input.bam', '--output-dir', str(tmp_path), '--variant-index', '2'])
    assert selected == [variant]
    assert exported[0]['inputs']['variant_index'] == 2
    assert exported[0]['inputs']['loaded_variant_count'] == 2


@pytest.mark.parametrize('index', ['0', '3'])
def test_plot_rejects_invalid_variant_index_before_opening_bam(monkeypatch, index):
    from isovar.cli import isovar_plot
    monkeypatch.setattr(isovar_plot, 'variant_collection_from_args', lambda args: [1, 2])
    monkeypatch.setattr(isovar_plot, 'alignment_file_from_args', lambda args: pytest.fail('BAM opened'))
    with pytest.raises(SystemExit) as error:
        isovar_plot.run(['--bam', 'input.bam', '--variant-index', index])
    assert error.value.code == 2


def test_evidence_report_uses_original_mates_from_producer_support():
    from types import SimpleNamespace
    from isovar import AlleleRead, ProteinSequenceCreator, ReadCollector
    from isovar.cli.isovar_evidence_report import collect

    variant = Variant('1', 1000, 'G', 'C', genome='GRCh38')
    placement = (0, 996, '7M', False)
    merged = AlleleRead('AAA', 'C', 'CCC', 'pair', source_read_count=2,
        source_alignments=((('rg', 'pair', 64), placement), (('rg', 'pair', 128), placement)))
    evidence = ReadEvidence(1000, 'G', 'C', [], [merged], [])
    class Collector(ReadCollector):
        def read_evidence_for_variant(self, variant, alignment):
            return evidence
    alignment = SimpleNamespace(header={'SQ': [{'SN': '1', 'LN': 2000}]})
    calls, _, _ = collect([variant], alignment, Collector(), ProteinSequenceCreator(),
                          sample_id='sample', source='library')
    assert calls[0]['alt_reads'] == calls[0]['total_reads'] == '2'
    assert calls[0]['alt_fragments'] == calls[0]['total_fragments'] == '1'
    assert calls[0]['ref_reads'] == calls[0]['other_reads'] == '0'


def test_protein_group_label_includes_covered_ranks(evidence_data, monkeypatch):
    import isovar.protein_comparison as comparison
    data = evidence_data[0]
    protein = data['modes'][0]['protein']
    rows = [dict(label='sample', source='sample', assembly=True, rank=2, covered_ranks=[1, 2],
                 protein=protein, frames=[('+', 0)], transcript_ids=['TX1']),
            dict(label='Varcode', source='Varcode', assembly=None, rank=1,
                 protein=protein, transcript_ids=['TX1'], predictions=[])]
    monkeypatch.setattr(comparison, 'comparison_rows', lambda products: rows)
    _, figure = next(comparison.protein_comparison_figures([dict(visualization=data)]))
    text = '\n'.join(t.get_text() for t in figure.axes[0].texts)
    assert 'ranks 1, 2' in text
