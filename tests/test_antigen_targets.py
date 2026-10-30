"""One target selection and reference preparation contract for all entry points."""

from contextlib import ExitStack
import json
from types import SimpleNamespace

import pandas as pd
import pysam
import pytest

from isovar import antigen_gene_selection, prepare_antigen_reference, reconstruct_antigen_isoforms
from isovar.antigen_targets import _annotation_genes
from isovar.cli.isovar_isoforms import run
from isovar.cohort import _antigen_path_task, _report_identity, _selected_file_keys
from tests.test_antigen_isoforms import reference_fixture


@pytest.fixture
def catalogue(monkeypatch):
    from oncoref import cta_selection

    monkeypatch.setattr(cta_selection, 'coverage', lambda **kw: dict(
        selection_status='available', policy_id='test-policy', reference={'reference_id': 'test-reference'}))
    monkeypatch.setattr(cta_selection, 'results', lambda **kw: pd.DataFrame([
        dict(gene_ids=['G', 'MISSING'], source_snapshot_id='sha256:test-snapshot',
             policy_id='test-policy', category='CTA', entity_id='test-group')]))
    monkeypatch.setattr(cta_selection, 'sources', lambda **kw: pd.DataFrame([
        dict(source_id='test-source', source_url='https://example.org/source', status='used')]))


def test_published_catalogue_retains_all_loci_policy_and_sources():
    selection = antigen_gene_selection(gene_set='oncoref-cta')
    assert selection['genes'] == sorted({g for row in selection['decisions'] for g in row['gene_ids']})
    assert selection['coverage']['policy_id']
    assert selection['coverage']['reference']['reference_id']
    assert all(row['source_snapshot_id'] for row in selection['decisions'])
    assert selection['provider']['version']
    assert selection['sources']
    json.dumps(selection, allow_nan=False)


@pytest.mark.parametrize('kwargs, error', [
    ({}, ValueError), ({'genes': [], 'gene_set': 'oncoref-cta'}, ValueError),
    ({'genes': []}, ValueError), ({'genes': 'PRAME'}, TypeError),
    ({'genes': ['']}, ValueError), ({'gene_set': 'unknown'}, ValueError),
])
def test_invalid_selection_is_rejected(kwargs, error):
    with pytest.raises(error):
        antigen_gene_selection(**kwargs)


def test_missing_catalogue_id_is_not_reinterpreted_as_symbol():
    class Annotation:
        def gene_by_id(self, query):
            raise ValueError(query)

        def genes_by_name(self, query):
            return [SimpleNamespace(id='OTHER_GENE')]

    with pytest.warns(UserWarning, match='absent'), pytest.raises(ValueError, match='No selected'):
        _annotation_genes(Annotation(), dict(kind='catalogue', genes=['MISSING_ID']))
    resolved, mapping = _annotation_genes(Annotation(), dict(kind='explicit', genes=['SYMBOL']))
    assert resolved == ['OTHER_GENE']
    assert mapping == [dict(query='SYMBOL', gene_ids=['OTHER_GENE'], status='gene_name')]


def test_api_cli_cohort_keep_same_models_and_selection(tmp_path, catalogue):
    fasta, gtf, bam_path = reference_fixture(tmp_path)
    reference = dict(genome_fasta=str(fasta), gtf=str(gtf), annotation_version='test',
                     alignment_reference='custom', reference_cache=str(tmp_path / 'cache'))
    with pysam.AlignmentFile(bam_path) as bam, pytest.warns(UserWarning, match='MISSING'):
        prepared = prepare_antigen_reference(bam, gene_set='oncoref-cta', **reference)
        expected = reconstruct_antigen_isoforms(
            bam, targets=prepared.targets, reference_identity=prepared.identity,
            sample_id='sample', source=str(bam_path), discovery_partition='gene')
    assert len(prepared.targets) == 2
    assert prepared.target_selection['gene_mapping'] == [
        dict(query='G', gene_ids=['G'], status='gene_id'),
        dict(query='MISSING', gene_ids=[], status='absent_from_annotation')]
    out = tmp_path / 'cli'
    with pytest.warns(UserWarning, match='MISSING'):
        run(['reconstruct', '--bam', str(bam_path), '--sample-id', 'sample',
             '--gene-set', 'oncoref-cta', '--genome-fasta', str(fasta), '--gtf', str(gtf),
             '--annotation-version', 'test', '--alignment-reference', 'custom',
             '--discovery-partition', 'gene',
             '--reference-cache', reference['reference_cache'], '--output-dir', str(out)])
    cli = json.loads((out / 'isoforms.json').read_text())
    files = dict(bam=bam_path, fasta=fasta, gtf=gtf)
    task = dict(gene_set='oncoref-cta', reference=dict(reference, genome_fasta='fasta', gtf='gtf'),
                parameters=dict(discovery_partition='gene'))
    source = dict(source_id=str(bam_path), sample_id='sample', bam='bam')
    with ExitStack() as contexts, pytest.warns(UserWarning, match='MISSING'):
        cohort = _antigen_path_task('isoforms', task, source, files, contexts, tmp_path, {})
    for result in (cli, cohort):
        assert result['target_selection'] == prepared.target_selection
        for key in ('isoforms', 'partial_sequences', 'protein_groups', 'support', 'reference', 'discovery'):
            assert json.loads(json.dumps(result[key])) == json.loads(json.dumps(expected[key]))


def test_cohort_reference_files_and_selection_are_part_of_receipt_identity():
    source = dict(source_id='rna', bam='bam')
    task = dict(source_id='rna', gene_set='oncoref-cta',
                reference=dict(genome_fasta='fasta', gtf='gtf', annotation_version='test'))
    manifest = dict(sources=[source], isoforms=[task],
                    files={key: dict(sha256=key) for key in ('bam', 'fasta', 'gtf')})
    assert _selected_file_keys(manifest, None) == {'bam', 'fasta', 'gtf'}
    first = _report_identity(task, source, manifest, 'isoforms')
    manifest['files']['gtf']['sha256'] = 'changed'
    assert first != _report_identity(task, source, manifest, 'isoforms')
    manifest['files']['gtf']['sha256'] = 'gtf'
    changed = dict(task, genes=['PRAME'])
    changed.pop('gene_set')
    assert first != _report_identity(changed, source, manifest, 'isoforms')


def test_named_catalogue_rejects_nonhuman_assembly(tmp_path, catalogue):
    with pytest.raises(ValueError, match='human CTA'):
        prepare_antigen_reference(None, gene_set='oncoref-cta', genome_fasta='unused', gtf='unused',
                                  annotation_version='test', assembly='GRCm39')
