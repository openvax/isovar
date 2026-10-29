"""Whole assessments must agree after one source-wide routing pass."""

from copy import deepcopy

import pysam
import pytest

from isovar.read_collector import ReadCollector
from isovar.rna_library import prepare_rna_library
from isovar.sv_path_routing import PathEvidenceQuery, _SeedRouter, route_path_evidence
from isovar.sv_rna_assessment import assess_sv_rna_paths
from isovar.sv_rna_path_support import reverse_complement
from tests.test_sv_rna import Scenario, record, write_bam
from tests.test_sv_rna_audit import manifest_fixture


@pytest.mark.parametrize('strand', ['+', '-'])
@pytest.mark.parametrize('minimum', [7, 20])
def test_routed_assessment_equals_exhaustive_with_all_original_alternatives(tmp_path, strand, minimum):
    scenario = Scenario(strand)
    originals = scenario.reads()
    originals += [record('f20', '1', 15000, '100M', scenario.sequence[200:300], flag=256),
                  record('f20', '1', 16000, '100M', scenario.sequence[200:300], flag=2048),
                  record('f20', '1', 15000, '100M', scenario.sequence[200:300], group='different'),
                  record('unplaced', '1', 0, None, scenario.sequence[190:290], flag=4),
                  record('remote', '2', 18000, '100M', scenario.sequence[180:280]),
                  record('no-quality', '2', 18000, '100M', scenario.sequence[180:280]),
                  record('raw-only', '2', 18000, '100M', 'A' * 100)]
    originals[-2].query_qualities = None
    originals += originals[:3]  # physical duplicates remain in source accounting
    path = write_bam(tmp_path / 'original.bam', originals)
    result = scenario.run(path)
    empty = dict(deepcopy(result), paths=[])
    other = deepcopy(result)
    other['paths'] = [dict(deepcopy(result['paths'][0]), path_id='second-branch')]
    results = {'first': result, 'empty': empty, 'other': other}
    references = (scenario.donor_ref, scenario.acceptor_ref)
    raw = [('a', 'raw-only', 0, scenario.sequence[170:300]),
           ('raw', 'unmapped-fastq', 2, reverse_complement(scenario.sequence[200:300]))]
    expected = {}
    with pysam.AlignmentFile(path) as bam:
        for key, value in results.items():
            expected[key] = assess_sv_rna_paths(value, bam, references, raw_reads=raw, min_path_match=minimum)
    # Index eligibility is deliberately narrower than path assessment.
    prepared = prepare_rna_library(path, tmp_path / 'cache',
        read_collector=ReadCollector(use_secondary_alignments=False))
    queries = tuple(PathEvidenceQuery(key, tuple(p['sequence'] for p in value['paths']), references, minimum)
                    for key, value in results.items())
    with prepared.open() as library, route_path_evidence(library, queries, tmp_path / 'routes.sqlite', raw_reads=iter(raw)) as routed:
        for key, value in results.items():
            actual = assess_sv_rna_paths(value, None, references, min_path_match=minimum, evidence=routed.evidence(key))
            assert actual == expected[key]


def test_pathless_queries_do_not_scan_sequences_and_preserve_all_exon_counts(tmp_path, monkeypatch):
    scenario = Scenario()
    path = write_bam(tmp_path / 'original.bam', scenario.reads())
    result = dict(scenario.run(path), paths=[])
    references = (scenario.donor_ref, scenario.acceptor_ref)
    with pysam.AlignmentFile(path) as bam:
        expected = assess_sv_rna_paths(result, bam, references)
    with prepare_rna_library(path, tmp_path / 'cache').open() as library:
        def no_scan():
            raise AssertionError('An empty path inventory must not scan original sequences')
        monkeypatch.setattr(library, 'original_records', no_scan)
        query = PathEvidenceQuery('empty', (), references)
        with route_path_evidence(library, (query,), tmp_path / 'routes.sqlite') as routed:
            actual = assess_sv_rna_paths(result, None, references, evidence=routed.evidence('empty'))
        assert actual == expected


@pytest.mark.parametrize('changed', ['sequence', 'match_threshold', 'exon_threshold', 'references'])
def test_reused_routed_evidence_rejects_changed_assessment_scope(tmp_path, changed):
    scenario = Scenario()
    path = write_bam(tmp_path / 'original.bam', scenario.reads())
    result = scenario.run(path)
    references = (scenario.donor_ref, scenario.acceptor_ref)
    query = PathEvidenceQuery('view', tuple(p['sequence'] for p in result['paths']), references)
    with prepare_rna_library(path, tmp_path / 'cache').open() as library:
        with route_path_evidence(library, (query,), tmp_path / 'routes.sqlite') as routed:
            evidence = routed.evidence('view')
    options = {}
    if changed == 'sequence':
        result = deepcopy(result)
        result['paths'][0]['sequence'] += 'A'
    elif changed == 'match_threshold':
        options['min_path_match'] = 19
    elif changed == 'exon_threshold':
        options['min_exon_overlap'] = 14
    else:
        references = references[:1]
    with pytest.raises(ValueError, match='evidence.*differ from assessment'):
        assess_sv_rna_paths(result, None, references, evidence=evidence, **options)


def test_assessment_requires_a_source_or_collected_evidence():
    with pytest.raises(ValueError, match='requires an original BAM'):
        assess_sv_rna_paths(dict(reference_name='GRCh38', paths=[]), None, ())


@pytest.mark.parametrize('sequence', ['ACGTN' * 60, 'NNNNNNN', 'acgtNNNNacgt', 'AA\0TT', 'éACTG'])
def test_seed_router_keeps_literal_sequence_semantics(sequence):
    paths = ('ACGTNACGTNNNNNNNACGT', 'AA\0TTéACTG')
    queries = tuple(PathEvidenceQuery(str(index), (path,), (), min_path_match=3) for index, path in enumerate(paths))
    router = _SeedRouter(queries)
    expected = {query.identifier for query in queries
                if any(word in sequence.upper() or reverse_complement(word) in sequence.upper()
                       for path in query.sequences for word in (path[i:i + 3] for i in range(len(path) - 2)))}
    assert router.matches(sequence) == expected


def test_fixed_width_numpy_router_equals_literal_oracle():
    import random
    randomizer = random.Random(71)
    paths = tuple(''.join(randomizer.choices('ACGTN', k=400)) for _ in range(4))
    queries = tuple(PathEvidenceQuery(str(index), (path,), (), min_path_match=9) for index, path in enumerate(paths))
    router = _SeedRouter(queries)
    for _ in range(80):
        sequence = ''.join(randomizer.choices('ACGTN', k=200))
        selected = randomizer.randrange(4)
        sequence += reverse_complement(paths[selected][30:70])
        expected = {str(index) for index, path in enumerate(paths)
                    if any(word in sequence or reverse_complement(word) in sequence
                           for word in (path[i:i + 9] for i in range(len(path) - 8)))}
        assert router.matches(sequence) == expected


@pytest.mark.parametrize('minimum', [1, 7, 15, 20, 31])
def test_spaced_seeds_cannot_drop_any_qualifying_exact_interval(minimum):
    import random
    randomizer = random.Random(912)
    path = ''.join(randomizer.choices('ACGTN', k=140))
    query = PathEvidenceQuery('path', (path,), (), min_path_match=minimum)
    router = _SeedRouter((query,))
    for start in range(len(path) - minimum + 1):
        window = path[start:start + minimum]
        assert router.matches('NN' + window + 'NN') == {'path'}
        assert router.matches('NN' + reverse_complement(window) + 'NN') == {'path'}


def multiple_view_manifest(tmp_path):
    from isovar.sv_rna_audit import write_json
    path, manifest = manifest_fixture(tmp_path)
    manifest['parameters'] = dict(dense_support=True, max_records=4, assemble=False)
    manifest['events'][0]['views'].append(dict(deepcopy(manifest['events'][0]['views'][0]), view_id='another-view'))
    manifest['events'].append(dict(deepcopy(manifest['events'][0]), event_id='another-event'))
    write_json(path, manifest)
    return path


def test_inventory_routes_one_original_pass_and_matches_exhaustive_assessment(tmp_path, monkeypatch):
    from contextlib import contextmanager
    from isovar import sv_rna_audit
    from isovar.rna_library import RnaLibraryIndex
    from isovar.sv_rna_assessment import collect_path_evidence
    manifest = multiple_view_manifest(tmp_path)
    scans = []
    original = RnaLibraryIndex.original_records

    def counted(library):
        scans.append(library.prepared.source)
        yield from original(library)

    monkeypatch.setattr(RnaLibraryIndex, 'original_records', counted)
    first = sv_rna_audit.run_sv_rna_audit(manifest, tmp_path / 'routed')
    assert len(scans) == 1

    @contextmanager
    def exhaustive(library, queries, database, raw_reads, read_collector):
        class Oracle:
            def evidence(self, identifier):
                query = next(query for query in queries if query.identifier == identifier)
                with pysam.AlignmentFile(library.prepared.source) as bam:
                    return collect_path_evidence(bam, query.references, (), query.min_exon_overlap,
                                                 read_collector=read_collector)
        yield Oracle()

    monkeypatch.setattr(sv_rna_audit, 'route_path_evidence', exhaustive)
    second = sv_rna_audit.run_sv_rna_audit(manifest, tmp_path / 'exhaustive')
    assert first == second
    for row in first['runs']:
        assert sv_rna_audit.read_json(tmp_path / 'routed' / row['artifact']) == sv_rna_audit.read_json(
            tmp_path / 'exhaustive' / row['artifact'])


def test_interrupted_inventory_resumes_verified_reconstructions(tmp_path, monkeypatch):
    from isovar import sv_rna_audit
    manifest = multiple_view_manifest(tmp_path)
    original = sv_rna_audit.route_path_evidence

    def interrupted(*args, **kwargs):
        raise OSError('Injected routing interruption')

    monkeypatch.setattr(sv_rna_audit, 'route_path_evidence', interrupted)
    output = tmp_path / 'output'
    with pytest.raises(OSError, match='Injected'):
        sv_rna_audit.run_sv_rna_audit(manifest, output)
    assert len(list((output / '.reconstructions').glob('*.json.gz'))) == 4
    monkeypatch.setattr(sv_rna_audit, 'route_path_evidence', original)

    def no_reconstruction(*args, **kwargs):
        raise AssertionError('Verified staged reconstruction should be reused')

    monkeypatch.setattr(sv_rna_audit, 'reconstruct_sv_rna', no_reconstruction)
    result = sv_rna_audit.run_sv_rna_audit(manifest, output, resume=True)
    assert len(result['runs']) == 8
    assert not list((output / '.reconstructions').iterdir())


def test_corrupted_reconstruction_checkpoint_is_not_accepted(tmp_path, monkeypatch):
    from isovar import sv_rna_audit
    manifest = multiple_view_manifest(tmp_path)

    def interrupted(*args, **kwargs):
        raise OSError('Injected routing interruption')

    monkeypatch.setattr(sv_rna_audit, 'route_path_evidence', interrupted)
    output = tmp_path / 'output'
    with pytest.raises(OSError):
        sv_rna_audit.run_sv_rna_audit(manifest, output)
    target = next((output / '.reconstructions').glob('*.json.gz'))
    target.write_bytes(b'corrupt')
    with pytest.raises(ValueError, match='Cannot resume changed'):
        sv_rna_audit.run_sv_rna_audit(manifest, output, resume=True)
