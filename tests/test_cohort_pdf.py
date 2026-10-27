"""Sequence display preserves codon anchors and shows deleted residues explicitly."""

from dataclasses import asdict, replace

import pytest

from isovar.cohort_pdf import aligned_columns, _fusion_comparisons
from isovar.fusion import reconstruct_fusion
from .test_fusion import example


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
    text = '\n'.join(page.extract_text() for page in PdfReader(report).pages)
    assert 'unavailable (no supported coding acceptor model)' in text
    assert 'MAAKGPFE' in (tmp_path / 'out/figures/sample/fusions/fusion/proteins.fasta').read_text()


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
