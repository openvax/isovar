"""One original-read policy must govern discovery, support and cached reuse."""

from contextlib import ExitStack
from dataclasses import asdict

import pysam
import pytest

from isovar import reconstruct_antigen_isoforms, reconstruct_viral_sequences
from isovar.cli.commands import run
from isovar.cohort import _antigen_path_task, run_cohort
from isovar.read_collector import ReadCollector
from isovar.read_end_inference import Adapter, ReadEndProfile
from isovar.rna_library import prepare_rna_library
from isovar.sv_rna_audit import file_sha256, read_json, write_json
from .test_antigen_isoforms import target_models, reference_fixture
from .test_sv_rna import aligned, record, write_bam


@pytest.fixture
def source(tmp_path):
    targets, positions = target_models()
    sequence = targets[0].reference.sequence
    reads = []
    for name, flag, mapq, qualities in [('good', 0, 60, True), ('low', 0, 0, True),
                                      ('duplicate', 1024, 60, True), ('secondary', 256, 60, True),
                                      ('no-quality', 0, 60, False)]:
        originals = aligned(name, sequence, positions[0], flag=flag, qualities=qualities)
        for read in originals:
            read.mapping_quality = mapq
        reads.extend(originals)
    path = write_bam(tmp_path / 'input.bam', reads)
    bases = ['N'] * 20000
    for (_, position, _), base in zip(positions[0], sequence):
        bases[position] = base
    fasta = tmp_path / 'genome.fa'
    fasta.write_text('>1\n' + ''.join(bases) + '\n>2\n' + 'N' * 20000 + '\n')
    return path, fasta, targets


def reconstruct(family, source, tmp_path, **kwargs):
    path, fasta, targets = source
    with pysam.AlignmentFile(path) as bam:
        common = dict(sample_id='sample', source='library', scratch_dir=tmp_path, **kwargs)
        if family == 'isoforms':
            return reconstruct_antigen_isoforms(bam, targets=targets, **common)
        return reconstruct_viral_sequences(bam, genome_fasta=fasta, contigs=['1'],
            reference_profile=None, reference_cache=tmp_path / 'reference-cache', **common)


@pytest.mark.parametrize('family', ['isoforms', 'viral'])
@pytest.mark.parametrize('settings,expected,excluded', [
    ({}, 4, {'duplicate'}),
    ({'min_mapping_quality': 20}, 3, {'duplicate', 'low'}),
    ({'use_duplicate_reads': True}, 5, set()),
    ({'use_secondary_alignments': False}, 3, {'duplicate', 'secondary'}),
    ({'use_reads_without_base_qualities': False}, 3, {'duplicate', 'no-quality'}),
])
def test_policy_controls_discovery_and_original_support(family, source, tmp_path, settings, expected, excluded):
    result = reconstruct(family, source, tmp_path, read_collector=ReadCollector.for_rna_paths(settings))
    assert result['support']['target_fragments'] == expected
    names = {row['fragment'][1] for row in result['support']['witnesses']}
    assert names == {'good', 'low', 'duplicate', 'secondary', 'no-quality'} - excluded
    discovery = result['discovery']
    selected = discovery['selected_segments'] if family == 'isoforms' else discovery['contigs'][0]['selected_segments']
    assert selected == expected
    assert all(result['collection'][name] == value for name, value in settings.items())


@pytest.mark.parametrize('family', ['isoforms', 'viral'])
def test_reused_index_matches_private_index_and_rejects_changed_policy(family, source, tmp_path):
    collector = ReadCollector.for_rna_paths()
    expected = reconstruct(family, source, tmp_path, read_collector=collector)
    reference = 'test' if family == 'isoforms' else None
    with prepare_rna_library(source[0], tmp_path / 'indexes', read_collector=collector,
                             reference_name=reference).open() as library:
        assert reconstruct(family, source, tmp_path, read_collector=collector, library_index=library) == expected
        with pytest.raises(ValueError, match='eligibility'):
            reconstruct(family, source, tmp_path, library_index=library,
                        read_collector=ReadCollector.for_rna_paths({'min_mapping_quality': 20}))
        with pytest.raises(ValueError, match='custom eligibility'):
            reconstruct(family, source, tmp_path, library_index=library,
                        read_collector=ReadCollector.for_rna_paths({'read_filter': lambda read: True}))


@pytest.mark.parametrize('family', ['isoforms', 'viral'])
def test_adapter_policy_changes_paths_without_changing_original_records(family, source, tmp_path):
    adapter = 'AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC'
    body = source[2][0].reference.sequence[40:140]
    path = write_bam(tmp_path / 'adapter.bam', [record('adapter', '1', 140, '100M34S', body + adapter)])
    source = path, source[1], source[2]
    collector = ReadCollector.for_rna_paths({'use_soft_clipped_bases': True})
    original = reconstruct(family, source, tmp_path, read_collector=collector)
    profile = ReadEndProfile('explicit-kit', adapters=(Adapter('adapter', adapter),))
    collector = ReadCollector.for_rna_paths({'use_soft_clipped_bases': True, 'read_end_profile': profile,
                                           'trim_adapters': True})
    trimmed = reconstruct(family, source, tmp_path, read_collector=collector)
    field = 'partial_sequences' if family == 'isoforms' else 'sequences'
    assert [row['sequence'] for row in original[field]] == [body + adapter]
    assert [row['sequence'] for row in trimmed[field]] == [body]
    assert original['support']['target_fragments'] == trimmed['support']['target_fragments'] == 1
    assert original['support']['witnesses'][0]['records'] == trimmed['support']['witnesses'][0]['records']
    assert trimmed['collection']['read_end_profile_sha256'] == profile.fingerprint


@pytest.mark.parametrize('family', ['isoforms', 'viral'])
def test_inapplicable_mate_merging_override_is_rejected(family, source, tmp_path):
    with pytest.raises(ValueError, match='merge_overlapping_fragments'):
        reconstruct(family, source, tmp_path,
                    read_collector=ReadCollector.for_rna_paths({'merge_overlapping_fragments': False}))


@pytest.mark.parametrize('family', ['isoforms', 'viral'])
def test_cli_uses_effective_policy_and_reuses_index(tmp_path, family):
    fasta, gtf, bam = reference_fixture(tmp_path)
    options = ['--bam', str(bam), '--genome-fasta', str(fasta), '--sample-id', 'sample',
               '--alignment-reference', 'custom', '--reference-cache', str(tmp_path / 'reference-cache'),
               '--library-cache', str(tmp_path / 'indexes'), '--min-mapping-quality', '61']
    command = ['viral', '--contig', 'chr1'] if family == 'viral' else [
        'isoforms', 'reconstruct', '--gene', 'CTA', '--gtf', str(gtf), '--annotation-version', 'test']
    filename = 'viral.json' if family == 'viral' else 'isoforms.json'
    for name in ['first', 'second']:
        run(command + options + ['--output-dir', str(tmp_path / name)])
        result = read_json(tmp_path / name / filename)
        assert result['support']['target_fragments'] == 0
        assert result['collection']['min_mapping_quality'] == 61
        assert result['input']['excluded_records'] == {'mapping_quality': 1}
    assert len(list((tmp_path / 'indexes').rglob('records.sqlite'))) == 1


@pytest.mark.parametrize('family', ['isoforms', 'viral'])
def test_cohort_applies_policy_and_pins_fasta_outputs(tmp_path, source, family):
    bam, fasta, targets = source
    target_path = tmp_path / 'targets.json'
    write_json(target_path, [asdict(t) for t in targets])
    task = dict(id='target', sample_id='sample', source_id='library',
                read_collection={'min_mapping_quality': 20})
    task.update(dict(targets='targets', reference_identity={}) if family == 'isoforms' else dict(
        genome_fasta='fasta', parameters=dict(contigs=['1'], reference_profile=None, reference_cache=str(tmp_path / 'cache'))))
    manifest = dict(schema='isovar.cohort_manifest.v1', samples=[dict(sample_id='sample')],
        sources=[dict(sample_id='sample', source_id='library', library_id='library', bam='bam',
                      scope='complete synthetic library', provenance={'fixture': 'original records'})],
        files={key: dict(path=str(path), sha256=file_sha256(path))
               for key, path in [('bam', bam), ('fasta', fasta), ('targets', target_path)]}, **{family: [task]})
    manifest_path = tmp_path / 'manifest.json'
    write_json(manifest_path, manifest)
    output = tmp_path / 'output'
    run_cohort(manifest_path, output, render=False)
    result = read_json(output / family / 'target/result.json.gz')
    assert result['support']['target_fragments'] == 3
    assert result['collection']['min_mapping_quality'] == 20
    receipt = read_json(output / family / 'target/receipt.json')
    assert 'proteins.fasta' in receipt['artifacts']
    run_cohort(manifest_path, output, resume=True, render=False)
    task['read_collection']['min_mapping_quality'] = 0
    write_json(manifest_path, manifest)
    with pytest.raises(ValueError, match='identity'):
        run_cohort(manifest_path, output, resume=True, render=False)


def test_cohort_reuses_library_across_views_and_separates_changed_policy(tmp_path, source, monkeypatch):
    import isovar.rna_library as module

    bam, fasta, targets = source
    target_path = tmp_path / 'targets.json'
    write_json(target_path, [asdict(t) for t in targets])
    prepared = []
    original = module.prepare_rna_library
    def tracked(*args, **kwargs):
        prepared.append(args[0])
        return original(*args, **kwargs)
    monkeypatch.setattr(module, 'prepare_rna_library', tracked)
    files = dict(bam=bam, fasta=fasta, targets=target_path)
    source_row = dict(source_id='library', sample_id='sample', bam='bam')
    task = dict(targets='targets', reference_identity={}, read_collection={'min_mapping_quality': 20})
    with ExitStack() as contexts:
        libraries = {}
        first = _antigen_path_task('isoforms', task, source_row, files, contexts, tmp_path, libraries)
        second = _antigen_path_task('isoforms', task, source_row, files, contexts, tmp_path, libraries)
        assert first == second
        assert len(prepared) == 1
        task['read_collection']['min_mapping_quality'] = 0
        third = _antigen_path_task('isoforms', task, source_row, files, contexts, tmp_path, libraries)
        assert third['support']['target_fragments'] == 4
        assert len(prepared) == 2
