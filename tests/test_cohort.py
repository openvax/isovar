"""The cohort runner reuses library reconstruction and preserves failed calls."""

import gzip
import json
from pathlib import Path

import pytest

from isovar.cohort import _literal, _reference_declaration, run_cohort, run_small_variant_source
from pyensembl import EnsemblRelease
from isovar.sv_rna_audit import file_sha256, write_json

from .testing_helpers import data_path


@pytest.mark.parametrize('source_contig', ['M', 'MT', 'chrM', 'chrMT'])
def test_mitochondrial_nomination_uses_annotation_alias_and_preserves_source(source_contig):
    variant = _literal(dict(contig=source_contig, start=12994, ref='G', alt='A'), EnsemblRelease(87))
    assert variant.contig == 'MT'
    assert variant.original_contig == source_contig
    assert variant.effects().top_priority_effect().gene_name == 'MT-ND5'


@pytest.mark.parametrize('source_contig,release,expected', [
    ('chr1', 87, '1'), ('chrUn_KI270442v1', 87, 'KI270442.1'),
    ('chrUn_gl000216', 75, 'GL000216.1'), ('GL000216.2', 75, 'GL000216.2'),
    ('chrUn_unrecognized', 87, 'chrUn_unrecognized')])
def test_contig_conversion_preserves_accession_versions_and_unknown_names(source_contig, release, expected):
    variant = _literal(dict(contig=source_contig, start=100, ref='A', alt='C'), EnsemblRelease(release))
    assert variant.contig == expected
    assert variant.original_contig == source_contig


def test_reference_declaration_recognizes_alias_but_rejects_different_build():
    genome = EnsemblRelease(87)
    matched = _reference_declaration({'reference_name': 'hg38'}, genome)
    assert matched['status'] == 'assembly_declaration_matches'
    assert matched['sequence_equivalence'] == 'not_established_by_assembly_name'
    with pytest.raises(ValueError, match='does not perform liftover'):
        _reference_declaration({'provenance': {'assembly': 'GRCh37'}}, genome)


def test_real_rna_resume_export_and_changed_inputs(tmp_path):
    bam = Path(data_path('data/colo829_whole_run/SRR8615617.all_sites.grch38.bam'))
    calls = tmp_path / 'calls.json'
    write_json(calls, [dict(contig='14', start=92811641, ref='GAGA', alt='G'),
                       dict(contig='not_a_contig', start=100, ref='A', alt='G')])
    task = dict(variants='calls', ensembl_release=87,
                protein_settings=dict(protein_sequence_length=33, protein_sequence_preference='support'))
    source = dict(sample_id='COLO829', source_id='original', bam='bam')
    files = dict(calls=calls, bam=bam)
    output = tmp_path / 'result'
    receipt = run_small_variant_source(task, source, files, output)
    assert receipt['status'] == 'completed_with_errors'
    assert receipt['errors'] == 0
    assert receipt['unresolved_chromosomes'] == 1
    with gzip.open(output / receipt['results'], 'rt') as stream:
        good, bad = [json.loads(line) for line in stream]
    assert good['gene_names'] == ['GOLGA5']
    proteins = good['protein_hypotheses']['events'][0]['protein_hypotheses']
    assert proteins and proteins[0]['amino_acids']
    comparison = proteins[0]['translations'][0]['reference_comparison']
    assert comparison['focal_amino_acid_change']
    assert bad['status'] == 'unresolved_chromosome'
    assert bad['contig_resolution']['status'] == 'missing'
    assert 'protein_hypotheses' not in bad
    prior = file_sha256(output / 'variants.jsonl.gz')
    again = run_small_variant_source(task, source, files, output, resume=True)
    assert again['identity_key'] == receipt['identity_key']
    assert prior == file_sha256(output / 'variants.jsonl.gz')
    with pytest.raises(ValueError, match='Checkpoint exists'):
        run_small_variant_source(task, source, files, output)
    write_json(calls, [dict(contig='14', start=92811641, ref='GAGA', alt='G')])
    with pytest.raises(ValueError, match='Resume inputs'):
        run_small_variant_source(task, source, files, output, resume=True)


def test_source_selection_does_not_require_other_library_inputs(tmp_path):
    earlier = tmp_path / 'assessment.json'
    write_json(earlier, {'scope': 'original assessment'})
    manifest = dict(schema='isovar.cohort_manifest.v1', samples=[{'sample_id': 'sample'}],
        files={'earlier': dict(path=str(earlier), sha256=file_sha256(earlier)),
               'pending-bam': dict(path='not-downloaded.bam', sha256='0' * 64),
               'calls': dict(path='not-downloaded.json', sha256='0' * 64)},
        sources=[dict(source_id=name, sample_id='sample', library_id=name,
                      scope='declared input scope', provenance={'accession': name},
                      **({'bam': 'pending-bam'} if name == 'pending' else {}))
                 for name in ['earlier', 'pending']],
        prior_assessments=[dict(id='prior', source_id='earlier', sample_id='sample', input='earlier')],
        small_variants=[dict(source_id='pending', variants='calls', ensembl_release=87)])
    path = tmp_path / 'manifest.json'
    write_json(path, manifest)
    assert run_cohort(path, tmp_path / 'out', source_ids=['earlier'], render=False) == []
    with pytest.raises(ValueError, match='Requested source'):
        run_cohort(path, tmp_path / 'out', source_ids=['unknown'], render=False)
    earlier.write_text('{}')
    with pytest.raises(ValueError, match='checksum mismatch'):
        run_cohort(path, tmp_path / 'out', source_ids=['earlier'], render=False)


def test_report_sample_id_cannot_escape_output_directory(tmp_path):
    path = tmp_path / 'manifest.json'
    write_json(path, dict(schema='isovar.cohort_manifest.v1', samples=[{'sample_id': '../outside'}],
                         sources=[], files={}))
    with pytest.raises(ValueError, match='Sample IDs'):
        run_cohort(path, tmp_path / 'out', render_only=True)


@pytest.mark.parametrize('family', ['sv_audits', 'fusions', 'viral', 'isoforms', 'prior_assessments'])
@pytest.mark.parametrize('render', [False, True])
def test_unavailable_additional_source_does_not_require_missing_inputs(tmp_path, family, render):
    manifest = dict(schema='isovar.cohort_manifest.v1', samples=[{'sample_id': 'sample'}],
        files={}, sources=[dict(source_id='missing', sample_id='sample', library_id='library',
                               scope='Original RNA', provenance={'accession': 'missing'},
                               unavailable_reason='Source server did not provide the file')])
    manifest[family] = [dict(id='task', source_id='missing', sample_id='sample', input='not-present')]
    path = tmp_path / 'manifest.json'
    write_json(path, manifest)
    paths = run_cohort(path, tmp_path / 'out', render=render)
    assert len(paths) == int(render)
    if render:
        assert paths[0].exists()
    entry = json.loads((tmp_path / 'out/cohort.json').read_text())['additional'][0]
    assert entry['status'] == 'unavailable'
    assert entry['reason'] == manifest['sources'][0]['unavailable_reason']
    manifest[family][0]['sample_id'] = 'another-sample'
    write_json(path, manifest)
    with pytest.raises(ValueError, match='different samples'):
        run_cohort(path, tmp_path / 'out', render=False)
