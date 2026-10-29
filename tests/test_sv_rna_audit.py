"""The unified workflow preserves source identity, uncertainty and original witnesses."""

from copy import deepcopy
from dataclasses import asdict
from pathlib import Path

import pysam
import pytest

from isovar.sv_rna_assessment import assess_sv_rna_paths, coding_intervals
from isovar.sv_rna_audit import file_sha256, read_json, run_sv_rna_audit, write_json
from isovar.sv_rna_reconciliation import hypotheses_from_varcode, reconcile_sv_rna_hypotheses
from tests.test_sv_rna import Scenario, record, write_bam


def assess(scenario, bam_path, result=None):
    result = scenario.run(bam_path) if result is None else result
    with pysam.AlignmentFile(str(bam_path)) as bam:
        return assess_sv_rna_paths(result, bam, [scenario.donor_ref, scenario.acceptor_ref], max_chain_reads=30)


@pytest.mark.parametrize('strand', ['+', '-'])
def test_all_paths_have_separate_junction_exon_frame_and_sequence_evidence(tmp_path, strand):
    s = Scenario(strand)
    bam = write_bam(tmp_path / 'rna.bam', s.reads())
    result = s.run(bam)
    original = deepcopy(result)
    # An unresolved branch must not acquire the coding branch's frame.
    unresolved = deepcopy(result['paths'][0])
    unresolved.update(path_id='unresolved', translations=[], frame_status='unresolved')
    result['paths'].append(unresolved)
    value = assess(s, bam, result)
    first, second = value['paths']
    assert first['proteins'] and first['proteins'][0]['event_linked']
    assert first['exon_links'] and first['exact_exon_anchors']
    assert first['junctions'][0]['exact_sequence_support']['complete_fragments']
    assert first['full_sequence_support']['overlap_chain']
    assert not first['full_sequence_support']['complete_fragments']
    assert second['frame_status'] == 'unresolved' and not second['proteins']
    assert result['paths'][0] == original['paths'][0]
    # The reverse-strand CDS starts in the last genomic exon.
    intervals = coding_intervals(s.acceptor_ref)
    assert sum(b - a for a, b in intervals) == s.acceptor_ref.cds_end - s.acceptor_ref.cds_start
    assert intervals[0][0] > intervals[-1][0] if strand == '-' else intervals[0][0] < intervals[-1][0]


def test_duplicate_processed_records_do_not_add_fragments(tmp_path):
    s = Scenario()
    reads = s.reads()
    bam = write_bam(tmp_path / 'rna.bam', reads)
    result = s.run(bam)
    one = assess(s, bam, result)
    copies = write_bam(tmp_path / 'copies.bam', reads + reads)
    two = assess(s, copies, result)
    assert two['records_scanned'] == 2 * one['records_scanned']
    assert two['fragments_scanned'] == one['fragments_scanned']
    assert two['paths'] == one['paths']


def test_original_alternative_placement_is_retained_as_competitor(tmp_path):
    s = Scenario()
    reads = s.reads()
    # f20 spans the junction; this competing secondary alignment places its
    # exact sequence elsewhere. It is an alternative placement, not a new read.
    reads.append(record('f20', '1', 15000, '100M', s.sequence[200:300], flag=256))
    bam = write_bam(tmp_path / 'rna.bam', reads)
    value = assess(s, bam)
    assert value['alternative_placements']
    reconciled = reconcile_sv_rna_hypotheses(value)
    placements = [r for r in reconciled['hypotheses'] if r['hypothesis']['kind'] == 'alternative_original_alignment']
    assert placements and all(r['status'] == 'indistinguishable_in_observed_RNA' for r in placements)
    assert any(r['hypothesis']['source']['alignment']['fragment_id'] == ['a', 'f20'] for r in placements)


def comparison_fixture():
    window = 'ACGTCAGTGCAT' + 'TGCA' + 'CGATGACTACGA'
    alternative = window[:12] + 'AAAAAA' + window[-12:]
    assessment = dict(event_id='event', sample_id='sample', source='library', reference_name='test',
                      reference_models=[], limitations=[], parameters=dict(junction_flank=12),
                      paths=[dict(path_id='observed', sequence=window, junctions=[dict(
                          window_sequence=window, relation='breakpoint_junction',
                          exact_sequence_support=dict(complete_fragments=[['rg', 'fragment']]))])])
    hypotheses = [dict(hypothesis_id='varcode', kind='varcode_prediction', source='Varcode', sequence=alternative),
                  dict(hypothesis_id='unresolved', kind='varcode_prediction', source='Varcode', sequence=None),
                  dict(hypothesis_id='partial', kind='ordinary_transcript', source='annotation', sequence=window[:15]),
                  dict(hypothesis_id='shared', kind='ordinary_transcript', source='annotation', sequence=window)]
    return assessment, hypotheses


def test_competition_requires_informative_original_intervals_and_keeps_all_hypotheses():
    assessment, hypotheses = comparison_fixture()
    result = reconcile_sv_rna_hypotheses(assessment, hypotheses)
    states = {r['hypothesis']['hypothesis_id']: r['status'] for r in result['hypotheses']}
    assert states == {'varcode': 'contradicted_by_informative_evidence', 'unresolved': 'unassessed',
                      'partial': 'insufficient_RNA', 'shared': 'supported', 'rna:observed': 'supported'}
    assert not result['translation_observed'] and not result['somatic_causation_proven']
    # With no discriminated alternative, matching reference and RNA remain tied.
    tied = reconcile_sv_rna_hypotheses(assessment, hypotheses[1:])
    assert tied['hypotheses'][-1]['status'] == 'indistinguishable_in_observed_RNA'
    # A reconstructed sequence with no original full-window witness is not evidence.
    assessment['paths'][0]['junctions'][0]['exact_sequence_support']['complete_fragments'] = []
    unobserved = reconcile_sv_rna_hypotheses(assessment, hypotheses)
    assert all(r['status'] in ('insufficient_RNA', 'unassessed') for r in unobserved['hypotheses'])
    assert not any(e['contradicted_hypothesis_ids'] for e in unobserved['evidence'])


def test_coexisting_conflicting_paths_are_not_ranked_or_dropped():
    assessment, hypotheses = comparison_fixture()
    other = deepcopy(assessment['paths'][0])
    other.update(path_id='alternative', sequence=hypotheses[0]['sequence'])
    other['junctions'][0].update(window_sequence=other['sequence'])
    other['junctions'][0]['exact_sequence_support']['complete_fragments'] = [['rg', 'other-fragment']]
    assessment['paths'].append(other)
    result = reconcile_sv_rna_hypotheses(assessment, hypotheses)
    assert all(r['status'] == 'mixed_support_and_contradiction' for r in result['hypotheses']
               if r['hypothesis']['hypothesis_id'] in ('varcode', 'shared', 'rna:observed', 'rna:alternative'))
    assert len(result['hypotheses']) == 6


def test_varcode_adapter_preserves_unresolved_and_materialized_candidates():
    from types import SimpleNamespace
    import varcode
    effect = SimpleNamespace(mutant_transcript=SimpleNamespace(cdna_sequence='ATGAAATAA'),
                             mutant_protein_sequence='MK', transcript=SimpleNamespace(id='tx'),
                             short_description='synthetic')
    unresolved = SimpleNamespace(mutant_transcript=None, mutant_protein_sequence=None,
                                 transcript=None, short_description='unresolved')
    candidates = [varcode.EffectCandidate(effect=effect, source='test'),
                  varcode.EffectCandidate(effect=unresolved, source='test')]
    values = hypotheses_from_varcode([SimpleNamespace(candidates=candidates)])
    assert {v['sequence'] for v in values} == {'ATGAAATAA', None}
    assert all(v['source']['varcode_version'] == varcode.__version__ for v in values)


def manifest_fixture(tmp_path):
    s = Scenario()
    bam = write_bam(tmp_path / 'rna.bam', s.reads())
    refs = tmp_path / 'references.json'
    write_json(refs, [asdict(s.donor_ref), asdict(s.acceptor_ref)])
    manifest = dict(schema='isovar.sv_rna_audit_manifest.v1',
                    files={k: dict(path=p.name, sha256=file_sha256(p)) for k, p in [('bam', bam), ('refs', refs)]},
                    sources=[dict(source_id='tumor', sample_id='sample', library_id='library', culture='culture-1',
                                  role='tumor', provenance=dict(caller='synthetic'), bam='bam'),
                             dict(source_id='normal', sample_id='normal', library_id='normal-library',
                                  culture='culture-2', role='normal', provenance=dict(caller='synthetic'),
                                  unavailable_reason='original normal reads unavailable')],
                    events=[dict(event_id='D--A', sources=['tumor', 'normal'], provenance=dict(caller='synthetic'),
                                 views=[dict(view_id='forward', references='refs', input=dict(
                                     reference_name='test', donor=asdict(s.donor), acceptor=asdict(s.acceptor), regions=[]))])])
    path = tmp_path / 'manifest.json'
    write_json(path, manifest)
    return path, manifest


def test_manifest_end_to_end_is_source_separated_and_resume_checks_outputs(tmp_path):
    manifest, _ = manifest_fixture(tmp_path)
    output = tmp_path / 'output'
    first = run_sv_rna_audit(manifest, output)
    assert len(first['runs']) == 2
    tumor, normal = first['runs']
    assert tumor['paths_with_event_linked_protein'] == 1
    assert normal['status'] == 'unassessed' and normal['retained_paths'] == 0
    snapshot = read_json(output / tumor['artifact'])
    assert snapshot['assessment']['source'] == 'tumor'
    assert snapshot['reconciliation']['hypotheses']
    assert first == run_sv_rna_audit(manifest, output, resume=True)
    with pytest.raises(ValueError, match='already exists'):
        run_sv_rna_audit(manifest, output)
    write_json(output / tumor['artifact'], {})
    with pytest.raises(ValueError, match='Cannot resume'):
        run_sv_rna_audit(manifest, output, resume=True)


def test_dense_inventory_prepares_one_library_for_all_events_and_views(tmp_path, monkeypatch):
    from isovar import rna_library, sv_rna_audit
    manifest_path, manifest = manifest_fixture(tmp_path)
    manifest['parameters'] = dict(dense_support=True, max_records=4, assemble=False)
    original_view = manifest['events'][0]['views'][0]
    manifest['events'][0]['views'].append(dict(deepcopy(original_view), view_id='second-budget', parameters=dict(max_records=6)))
    manifest['events'].append(dict(deepcopy(manifest['events'][0]), event_id='another-nomination'))
    write_json(manifest_path, manifest)
    prepared = []
    original = sv_rna_audit.prepare_rna_library

    def counted(*args, **kwargs):
        value = original(*args, **kwargs)
        prepared.append(value.directory)
        return value

    monkeypatch.setattr(sv_rna_audit, 'prepare_rna_library', counted)
    cache = tmp_path / 'library-cache'
    output = tmp_path / 'first'
    first = run_sv_rna_audit(manifest_path, output, library_cache_dir=cache)
    assert len(prepared) == 1
    assert len(first['runs']) == 8
    assert sum(row['status'] == 'unassessed' for row in first['runs']) == 4

    def no_rebuild(*args, **kwargs):
        raise AssertionError('An unchanged library must reuse its verified index')

    monkeypatch.setattr(rna_library, '_index_records', no_rebuild)
    second_output = tmp_path / 'second'
    second = run_sv_rna_audit(manifest_path, second_output, library_cache_dir=cache)
    assert first == second
    for row in first['runs']:
        assert read_json(output / row['artifact']) == read_json(second_output / row['artifact'])
    assert len(prepared) == 2 and prepared[0] == prepared[1]
    assert first == run_sv_rna_audit(manifest_path, output, library_cache_dir=cache, resume=True)
    assert len(prepared) == 2  # Resume needs no reconstruction index at all.


def test_manifest_verifies_inputs_before_creating_output(tmp_path):
    path, manifest = manifest_fixture(tmp_path)
    manifest['files']['bam']['sha256'] = '0' * 64
    write_json(path, manifest)
    output = tmp_path / 'output'
    with pytest.raises(ValueError, match='checksum mismatch'):
        run_sv_rna_audit(path, output)
    assert not output.exists()


def test_supplied_fusion_uses_existing_validator_and_partial_candidate(tmp_path):
    from tests.testing_helpers import fusion_input
    manifest_path, manifest = manifest_fixture(tmp_path)
    from isovar.sid_data import restore_records
    data = fusion_input(restore_records(read_json(Path(__file__).parent / 'data/fusions/coding-corpus/BCR--ABL1.input.json.gz')))
    supplied = tmp_path / 'fusion.json'
    write_json(supplied, data)
    refs = tmp_path / 'supplied_refs.json'
    write_json(refs, data['references'])
    for key, path in [('supplied', supplied), ('supplied_refs', refs)]:
        manifest['files'][key] = dict(path=path.name, sha256=file_sha256(path))
    fusion = data['fusion']
    manifest['sources'][0]['sample_id'] = fusion['provenance']['sample_id']
    event = manifest['events'][0]
    event.update(event_id=fusion['event_id'], sources=['tumor'])
    view = event['views'][0]
    view.update(supplied_fusion='supplied', references='supplied_refs', input=dict(
        reference_name=fusion['reference_name'], donor=fusion['donor'], acceptor=fusion['acceptor'], regions=[]))
    write_json(manifest_path, manifest)
    output = tmp_path / 'out'
    report = run_sv_rna_audit(manifest_path, output)
    snapshot = read_json(output / report['runs'][0]['artifact'])
    assert snapshot['reconstruction']['schema'].startswith('isovar.fusion_rna.')
    assert snapshot['assessment']['paths'][0]['proteins']
    assert all(p['event_linked'] for p in snapshot['assessment']['paths'][0]['proteins'])


def test_reconciliation_is_order_independent_and_keeps_repeated_flanks_unresolved():
    assessment, hypotheses = comparison_fixture()
    window = assessment['paths'][0]['sequence']
    hypotheses[0]['sequence'] = window[:12] + hypotheses[0]['sequence']
    first = reconcile_sv_rna_hypotheses(assessment, hypotheses)
    second = reconcile_sv_rna_hypotheses(assessment, hypotheses[::-1])
    keyed = lambda result: {r['hypothesis']['hypothesis_id']: r for r in result['hypotheses']}
    assert keyed(first) == keyed(second)
    assert keyed(first)['varcode']['status'] == 'insufficient_RNA'


def test_manifest_rejects_ambiguous_identities_and_unknown_selections(tmp_path):
    path, manifest = manifest_fixture(tmp_path)
    with pytest.raises(ValueError, match='not in manifest'):
        run_sv_rna_audit(path, tmp_path / 'out', event_ids=['missing'])
    manifest['sources'].append(manifest['sources'][0])
    write_json(path, manifest)
    with pytest.raises(ValueError, match='Duplicate source'):
        run_sv_rna_audit(path, tmp_path / 'out')
    assert not (tmp_path / 'out').exists()


def test_original_read_hcc1395_complete_protein_through_shared_assessment():
    from reports.hcc1395_colo829.build_hcc1395_eif3k_fusion import report, load_references
    data = Path(__file__).parent / 'data/hcc1395_eif3k_cyp39a1'
    value = report(data / 'original_unmapped.bam', data / 'models.json.gz',
                   data / 'ont_candidate_sequences.fasta.gz', data / 'coding_tiling.bam',
                   data / 'ont_candidate_alignments.bam', data / 'ont_exon_models.json',
                   data / 'dna_tumor_loci.bam', data / 'dna_normal_loci.bam')
    from isovar.read_identity import read_group
    with pysam.AlignmentFile(str(data / 'original_unmapped.bam')) as originals:
        raw = [(read_group(r), r.query_name, 1 if r.is_read1 else 2 if r.is_read2 else 0, r.query_sequence)
               for r in originals.fetch(until_eof=True) if r.query_sequence]
    with pysam.AlignmentFile(str(data / 'coding_tiling.bam')) as bam:
        assessment = assess_sv_rna_paths(value['result'], bam, load_references(data / 'models.json.gz'),
                                        raw_reads=raw, min_overlap=20, max_chain_reads=30)
    path, = assessment['paths']
    candidate, = path['proteins']
    assert candidate['completeness'] == 'complete' and candidate['event_linked']
    assert len(candidate['candidate']['amino_acids']) == 531
    assert candidate['sequence_support']['overlap_chain']
    assert not assessment['translation_observed']
    reconciliation = reconcile_sv_rna_hypotheses(assessment)
    assert reconciliation['hypotheses'] and reconciliation['sample_id'] == 'HCC1395'


def test_original_fastq_names_are_explicit_and_duplicate_segments_rejected(tmp_path):
    from isovar.sv_rna_audit import raw_fastq_reads
    path = tmp_path / 'original.fq'
    path.write_text('@SRR1.5.1\nACGT\n+\nIIII\n@SRR1.5.2\nTGCA\n+\nIIII\n')
    assert list(raw_fastq_reads(path, read_group='rg')) == [('rg', 'SRR1.5', 1, 'ACGT'), ('rg', 'SRR1.5', 2, 'TGCA')]
    path.write_text(path.read_text() * 2)
    with pytest.raises(ValueError, match='Duplicate'):
        list(raw_fastq_reads(path))


def test_installed_cli_contract_does_not_return_report_as_exit_status(tmp_path):
    from isovar.cli.commands import run
    manifest, _ = manifest_fixture(tmp_path)
    assert run(['sv-audit', '--manifest', str(manifest), '--output-dir', str(tmp_path / 'out')]) is None
    assert (tmp_path / 'out/report.json').exists()


def test_predicted_transcript_loss_is_retained_without_treating_no_rna_as_validation():
    assessment, hypotheses = comparison_fixture()
    hypotheses.append(dict(hypothesis_id='exon-loss', kind='varcode_prediction', source='Varcode', sequence=''))
    result = reconcile_sv_rna_hypotheses(assessment, hypotheses)
    row = next(r for r in result['hypotheses'] if r['hypothesis']['hypothesis_id'] == 'exon-loss')
    assert row['status'] == 'unassessed' and row['sequence_state'] == 'predicted_empty_transcript'
    assert row['hypothesis']['sequence'] == ''
    assert not any(r['status'] == 'contradicted_by_informative_evidence' for r in row['comparisons'])


def test_committed_cohort_manifest_preserves_every_original_inventory_row():
    from reports.hcc1395_colo829.build_unified_sv_manifest import OUTPUT, DATA
    manifest = read_json(OUTPUT / 'manifest.json')
    source_rows = {line for name in ('colo829-sv-grch38.vcf', 'colo829-sv-2022-hg38lifted.vcf')
                   for line in (DATA / name).read_text().splitlines() if not line.startswith('#')}
    assert len(manifest['events']) == len(source_rows) == 188
    assert {e['provenance']['original_vcf_row'] for e in manifest['events']} == source_rows
    assert sum(bool(e.get('unassessed_reason')) for e in manifest['events']) == 14
    assert sum(len(e.get('views', [])) * len(e['sources']) if e.get('views') else len(e['sources'])
               for e in manifest['events']) == 380
