"""Sequence display preserves codon anchors and shows deleted residues explicitly."""

from dataclasses import asdict, replace
from collections import Counter
import gzip
import json

import pytest

from isovar.cohort_pdf import aligned_columns, _fusion_comparisons
from isovar.fusion import reconstruct_fusion
from .test_fusion import example


@pytest.mark.parametrize('depths', [[109, 221], [100, 101], [1, 2, 200, 201], [1, 10, 100, 1000000], [221]])
def test_log_depth_tick_labels_do_not_overlap(tmp_path, monkeypatch, depths):
    from matplotlib.figure import Figure
    from matplotlib.backends.backend_agg import FigureCanvasAgg
    from isovar.cohort_pdf import _category_landscape
    figures = []
    original = Figure.savefig

    def capture(figure, *args, **kwargs):
        figures.append(figure)
        return original(figure, *args, **kwargs)

    monkeypatch.setattr(Figure, 'savefig', capture)
    _category_landscape(dict(source_id='synthetic', pass_id='ticks'), {'SNV': len(depths)},
                        {'SNV': depths}, tmp_path)
    figure = figures[0]
    canvas = FigureCanvasAgg(figure)
    canvas.draw()
    axis = figure.axes[1]
    boxes = [label.get_window_extent(canvas.get_renderer()) for label in axis.get_xticklabels() if label.get_visible()]
    assert 1 <= len(boxes) <= 3
    assert all(left.x1 < right.x0 for left, right in zip(boxes, boxes[1:]))
    assert list(axis.get_xticks())[0] == min(depths)
    assert list(axis.get_xticks())[-1] == max(depths)
    assert not axis.xaxis.get_minorticklocs().size


def test_funnel_legend_does_not_cover_the_data(tmp_path, monkeypatch):
    from matplotlib.figure import Figure
    from matplotlib.backends.backend_agg import FigureCanvasAgg
    from isovar.cohort_pdf import _landscape
    figures = []
    original = Figure.savefig

    def capture(figure, *args, **kwargs):
        figures.append(figure)
        return original(figure, *args, **kwargs)

    monkeypatch.setattr(Figure, 'savefig', capture)
    _landscape('sample', [(dict(source_id='source', pass_id='default'),
                          Counter(nominated=2, alt=2, changed_protein=1), {}, {}, {})], tmp_path)
    canvas = FigureCanvasAgg(figures[0])
    canvas.draw()
    axes = figures[0].axes[0]
    legend = axes.get_legend().get_window_extent(canvas.get_renderer())
    assert legend.y1 < axes.get_window_extent().y0


@pytest.mark.parametrize('scopes,expected_count,expected_label', [
    ([{}, {'status': 'not_assessed'}], 0, 'Not assessed'),
    ([{'status': 'assessed', 'overlapping_models': 0}], 1, '0'),
    ([{'status': 'no_alt_reads', 'overlapping_models': 0}], 1, '0'),
    ([{'status': 'assessed', 'overlapping_models': 0}, {'status': 'not_assessed'}], 1, '0 (partial)'),
])
def test_ncorf_scope_keeps_measured_zero_separate_from_missing_catalogue(
        tmp_path, scopes, expected_count, expected_label):
    from isovar.cohort_pdf import _collect, _ncorf_count, _ncorf_scope_description
    from isovar.sv_rna_audit import file_sha256
    path = tmp_path / 'results.jsonl.gz'
    with gzip.open(path, 'wt') as stream:
        for index, scope in enumerate(scopes):
            row = dict(status='assessed', category='SNV', exon_overlap=False,
                       nomination=dict(contig='1', start=index+1, ref='A', alt='G'), nc_orfs=scope,
                       protein_hypotheses=dict(events=[dict(allele_support={'alt': {'fragments': 0}},
                                                            protein_hypotheses=[])]))
            stream.write(json.dumps(row) + '\n')
    result = dict(status='completed', output=str(tmp_path), results=path.name,
                  result_sha256=file_sha256(path), nominated_alleles=len(scopes))
    counts, _, _, _ = _collect(result)
    assert counts['nc_catalogue_assessed'] == expected_count
    assert _ncorf_count(counts, 'changed_nc') == expected_label
    assert _ncorf_count(counts, 'changed_nc_span20') == expected_label
    text = _ncorf_scope_description(counts)
    if not expected_count:
        assert 'unassessed' in text and '0 overlap' not in text
    else:
        assert f'{expected_count}/{len(scopes)} nominations' in text


def test_ncorf_unknown_scope_does_not_erase_positive_candidates():
    from isovar.cohort_pdf import _ncorf_count
    counts = Counter(nominated=3, changed_nc=1)
    assert _ncorf_count(counts, 'changed_nc') == '1 (scope unknown)'


@pytest.mark.parametrize('depths,expected_x,expected_y', [
    ([20], [20], [1]),
    ([5, 20, 5, 0], [5, 20], [3, 1]),
    ([8, 8, 8], [8], [3]),
    ([0, 0], [], []),
])
def test_support_landscape_has_visible_inclusive_threshold_counts(tmp_path, monkeypatch, depths, expected_x, expected_y):
    from matplotlib.figure import Figure
    from isovar.cohort_pdf import _category_landscape
    figures = []
    save = Figure.savefig

    def capture(figure, *args, **kwargs):
        figures.append(figure)
        return save(figure, *args, **kwargs)

    monkeypatch.setattr(Figure, 'savefig', capture)
    path = _category_landscape(dict(source_id='synthetic', pass_id='thresholds'),
                               {'SNV': len(depths)}, {'SNV': depths}, tmp_path)
    assert path.is_file()
    support = figures[0].axes[1]
    if expected_x:
        line, = support.lines
        assert list(line.get_xdata()) == expected_x
        assert list(line.get_ydata()) == expected_y
        assert line.get_marker() not in (None, '', 'None')
        assert line.get_drawstyle() == 'steps-pre'
    else:
        assert not support.lines
        assert any(text.get_text() == 'No positive RNA ALT support' for text in support.texts)


def test_deletion_in_repeat_keeps_declared_anchor_and_neutral_flanks():
    ref = 'LAERKYSDEKKRVDELQQQVKLYKLNLESSKQE'
    rna = 'LAERKYSDEKRVDELQQQVKLYKLNLESSKQE'
    columns = aligned_columns(ref, rna, [9, 9], allele_length_delta=3)
    deleted = [i for i, c in enumerate(columns) if c['focal_deletion']]
    assert deleted == [9]
    assert columns[9]['reference'] == 'K'
    assert columns[9]['protein'] == '-'
    assert all(c['reference'] == c['protein'] for c in columns if not c['focal_deletion'])
    assert ''.join(c['protein'] for c in columns if c['position'] is not None) == rna


def test_frameshift_is_not_aligned_back_to_a_spurious_matching_tail():
    columns = aligned_columns('MKTLLAG', 'MKAGLLA', [2, 7], allele_length_delta=1, frameshift=True)
    assert [c['position'] for c in columns] == list(range(7))
    assert [c['reference_position'] for c in columns] == list(range(7))
    assert not any(c['focal_deletion'] for c in columns)


def test_truncated_reference_does_not_invent_matching_stop():
    columns = aligned_columns('MKAAALLL', 'MK', [2, 2], stop=True, reference_stop=True)
    assert columns[-1]['protein'] == '*'
    assert columns[-1]['reference'] != '*'


@pytest.mark.parametrize('donor_strand', ['+', '-'])
@pytest.mark.parametrize('acceptor_strand', ['+', '-'])
@pytest.mark.parametrize('partial', [0, 1, 2, 3, 4, 5])
@pytest.mark.parametrize('insert', ['', 'GCT', 'GCTGCT'])
def test_fusion_parents_follow_actual_codon_and_insert_offsets(donor_strand, acceptor_strand, partial, insert):
    fusion, refs, reads = example(donor_strand, acceptor_strand, partial=partial, insert=insert)
    result = reconstruct_fusion(fusion, refs, reads)
    panels, fasta = _fusion_comparisons(result, {'references': [asdict(r) for r in refs]},
                                      [2] * len(fusion.sequence), 2)
    panel, = panels
    protein = panel['translation']['amino_acids']
    donor, acceptor = map(''.join, panel['parents'])
    donor_length = (fusion.junction_start - panel['translation']['translation_start']) // 3
    inserted_length = len(insert) // 3
    assert donor[:donor_length] == protein[:donor_length]
    assert donor[donor_length:].strip() == ''
    assert acceptor[:donor_length + inserted_length].strip() == ''
    assert acceptor[donor_length + inserted_length:] == 'GPFE'
    assert protein[donor_length + inserted_length:] == 'GPFE'
    assert panel['junction'] == [donor_length * 3, (donor_length + inserted_length) * 3]
    assert set(panel['depth']) == {2}
    assert fasta['fusion_1'] == protein


@pytest.mark.parametrize('acceptor_kind', ['missing', 'noncoding', 'utr'])
def test_fusion_without_coding_acceptor_keeps_observed_protein(acceptor_kind):
    fusion, refs, reads = example()
    if acceptor_kind == 'missing':
        refs = refs[:1]
    elif acceptor_kind == 'noncoding':
        refs = (refs[0], replace(refs[1], cds_start=None, cds_end=None))
    else:
        acceptor = replace(refs[1], sequence='ATGGGTATGTTTGAATAA', cds_start=6)
        refs = (refs[0], acceptor)
        fusion = replace(fusion, sequence=fusion.sequence[:12] + acceptor.sequence[3:])
        reads = tuple(replace(read, sequence=fusion.sequence) for read in reads)
    result = reconstruct_fusion(fusion, refs, reads)
    assert result['status'] == 'translated'
    panels, fasta = _fusion_comparisons(result, {'references': [asdict(r) for r in refs]},
                                      [2] * len(fusion.sequence), 2)
    panel, = panels
    expected = 'MAAKGMFE' if acceptor_kind == 'utr' else 'MAAKGPFE'
    assert panel['translation']['amino_acids'] == fasta['fusion_1'] == expected
    assert panel['acceptor'] is None
    assert ''.join(panel['parents'][0]) == 'MAAK    '
    assert ''.join(panel['parents'][1]).strip() == ''
    assert 'reference_acceptor.1' not in fasta


@pytest.mark.parametrize('acceptor_strand', ['+', '-'])
@pytest.mark.parametrize('cut,insert', [(13, ''), (14, ''), (12, 'A'), (12, 'TG')])
def test_split_acceptor_codon_uses_contributing_reference_bases(acceptor_strand, cut, insert):
    fusion, refs, reads = example(acceptor_strand=acceptor_strand, cut=cut, insert=insert)
    result = reconstruct_fusion(fusion, refs, reads)
    panels, _ = _fusion_comparisons(result, {'references': [asdict(r) for r in refs]},
                                   [2] * len(fusion.sequence), 2)
    panel, = panels
    first_acceptor_codon = (fusion.junction_end - panel['translation']['translation_start']) // 3
    assert panel['parents'][1][first_acceptor_codon] == 'G'
    assert panel['acceptor_in_frame'] is False


def test_donor_only_fusion_renders_complete_cohort_pdf(tmp_path):
    from pypdf import PdfReader
    from isovar.cohort import run_cohort
    from isovar.sv_rna_audit import file_sha256, read_json, write_json
    from .test_cohort_integrity import fusion_manifest

    manifest = fusion_manifest(tmp_path)
    supplied_path = tmp_path / 'fusion.json'
    supplied = read_json(supplied_path)
    supplied['references'] = supplied['references'][:1]
    write_json(supplied_path, supplied)
    manifest['files']['fusion']['sha256'] = file_sha256(supplied_path)
    path = tmp_path / 'manifest.json'
    write_json(path, manifest)
    report, = run_cohort(path, tmp_path / 'out')
    pages = [page.extract_text() for page in PdfReader(report).pages]
    text = '\n'.join(pages)
    assert 'unavailable (no supported coding acceptor model)' in text
    assert 'MAAKGPFE' in (tmp_path / 'out/figures/sample/fusions/fusion/proteins.fasta').read_text()
    provenance_page, = [page for page in pages if 'Source provenance' in page]
    assert 'synthetic regression' in provenance_page


def test_fusion_reference_comparison_uses_annotation_and_transcript_identity():
    from isovar.fusion import FusionReference

    fusion, refs, reads = example(partial=3, cut=9)
    alternative = FusionReference(refs[0].transcript_id, 'test', 'another annotation', 'd', '+',
        ((0, 4), (103, 109), (300, 305)), 'ATGA' + refs[0].sequence[3:9] + 'CCTAA', 0, 15)
    refs += (alternative,)
    result = reconstruct_fusion(fusion, refs, reads)
    panels, fasta = _fusion_comparisons(result, {'references': [asdict(r) for r in refs]},
                                      [2] * len(fusion.sequence), 2)
    first = next(p for p in panels if p['donor_annotation'] == 'test annotation 1')
    assert first['translation']['amino_acids'] == 'AAGPFE'
    assert ''.join(first['parents'][0]) == 'AA    '
    assert {v for k, v in fasta.items() if k.startswith('reference_donor.1')} == {'MAAKPE', 'MSCS'}


def test_ambiguous_acceptor_reference_is_not_selected_by_input_order():
    fusion, refs, reads = example()
    refs += (replace(refs[1], annotation='another annotation'),)
    for ordered in (refs, refs[::-1]):
        result = reconstruct_fusion(fusion, ordered, reads)
        panels, _ = _fusion_comparisons(result, {'references': [asdict(r) for r in ordered]},
                                       [2] * len(fusion.sequence), 2)
        assert panels[0]['acceptor_status'] == 'ambiguous'
        assert panels[0]['acceptor'] is None
        assert ''.join(panels[0]['parents'][1]).strip() == ''


def test_codon_boundary_without_acceptor_does_not_invent_mutant_residues():
    from reportlab.lib.colors import HexColor
    from isovar.cohort_pdf import FusionSequencePanel, INK, UNKNOWN, MUTANT

    class Canvas:
        def __init__(self):
            self.residues = []

        def setFillColor(self, color):
            self.color = color

        def drawString(self, x, y, text):
            if y == 45 and len(text) == 1:
                self.residues.append((text, self.color))

        def __getattr__(self, name):
            return lambda *args, **kwargs: None

    canvas = Canvas()
    panel = FusionSequencePanel('MAAKGPFE', ('MAAK    ', '        '), [2] * 8, [12, 12], 0)
    panel.canv = canvas
    panel.draw()
    assert canvas.residues[:4] == [(aa, HexColor(INK)) for aa in 'MAAK']
    assert canvas.residues[4:] == [(aa, HexColor(UNKNOWN)) for aa in 'GPFE']
    canvas = Canvas()
    shifted = FusionSequencePanel('MAAKGDA', ('MAAK   ', '    GPF'), [2] * 7, [12, 12], 0,
                                  acceptor_in_frame=False)
    shifted.canv = canvas
    shifted.draw()
    assert canvas.residues[4] == ('G', HexColor(INK))
    assert canvas.residues[5:] == [(aa, HexColor(MUTANT)) for aa in 'DA']


def test_same_event_in_two_libraries_keeps_distinct_figures_and_fastas(tmp_path, monkeypatch):
    from copy import deepcopy
    import isovar.cohort_pdf as reporting
    from isovar.cohort import run_cohort
    from isovar.sv_rna_audit import file_sha256, read_json, write_json
    from .test_cohort_integrity import fusion_manifest

    manifest = fusion_manifest(tmp_path)
    supplied = read_json(tmp_path / 'fusion.json')
    supplied['fusion']['event_id'] = 'donor / acceptor'
    write_json(tmp_path / 'fusion.json', supplied)
    manifest['files']['fusion']['sha256'] = file_sha256(tmp_path / 'fusion.json')
    second = deepcopy(supplied)
    second['reads'].append(dict(second['reads'][0], read_id='third', fragment_id='third'))
    for read in second['reads']:
        read['library_id'] = 'second'
    write_json(tmp_path / 'second.json', second)
    manifest['files']['second'] = dict(path=str(tmp_path / 'second.json'), sha256=file_sha256(tmp_path / 'second.json'))
    manifest['sources'].append(dict(manifest['sources'][0], source_id='second', library_id='second'))
    manifest['fusions'].append(dict(manifest['fusions'][0], id='second', source_id='second', input='second'))
    path = tmp_path / 'manifest.json'
    write_json(path, manifest)
    images = []
    original_image = reporting.Image

    def image(path, **kwargs):
        if str(path).endswith('coverage.png'):
            images.append(str(path))
        return original_image(path, **kwargs)

    monkeypatch.setattr(reporting, 'Image', image)
    run_cohort(path, tmp_path / 'out')
    assert len(images) == len(set(images)) == 2
    directories = [tmp_path / 'out/figures/sample/fusions' / task for task in ('fusion', 'second')]
    assert all((directory / 'proteins.fasta').exists() for directory in directories)
    assert file_sha256(directories[0] / 'coverage.png') != file_sha256(directories[1] / 'coverage.png')
