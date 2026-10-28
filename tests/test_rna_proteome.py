"""RNA-source adapters retain uncertain outcomes and source-scoped evidence."""

from copy import deepcopy
from dataclasses import replace
import gzip
from hashlib import sha256
import json

from orfdata import Proteome, ReferenceContext
import pytest

from isovar import export_sv_rna_orfs, from_rna_result, reconstruct_fusion
from isovar.cli import commands
from tests.test_fusion import example
from tests.test_proteome_adapter import exported, records_of, REFERENCE
from tests.test_sv_rna_orf_export import reconstruction
from tests.test_variant_rna_orfs import assess, noncoding_rna as noncoding_rna


SV_REFERENCE = ReferenceContext(9606, 'human', 'synthetic', 'test', '1')
FUSION_REFERENCE = ReferenceContext(9606, 'human', 'test', 'test annotation', '1')


def test_sv_native_and_offline_keep_all_placements_and_zero_support(tmp_path):
    result = reconstruction()
    other = deepcopy(result['paths'][0])
    other['path_id'] = 'unsupported'
    other['exploratory_orfs']['candidates'][0]['full_interval_support']['witnesses'] = []
    result['paths'].append(other)
    before = deepcopy(result)
    export = export_sv_rna_orfs(result)
    native = from_rna_result(result, reference=SV_REFERENCE)
    offline = from_rna_result(json.loads(json.dumps(result)), reference=SV_REFERENCE)
    assert native.identity == offline.identity and result == before
    compact = from_rna_result(export, reference=SV_REFERENCE)
    assert [p.accession for p in native] == [p.accession for p in compact]
    protein, = native
    assert protein.sequence == 'MKMK'
    assert protein.producer.origin == 'sv_orf'
    assert protein.producer.start_complete is None and protein.producer.stop_complete is True
    assert not protein.evidence and protein.producer.tool_version is None
    placements = records_of(native, 'occurrence')
    assert {p.payload['path_id']: p.payload['rna_support']['fragments'] for p in placements} == {
        'path1': 1, 'unsupported': 0}
    hypothesis, = records_of(native, 'hypothesis')
    assert hypothesis.status == 'exploratory'
    assert hypothesis.payload['translation_observed'] is False
    assert hypothesis.payload['initiation_observed'] is False
    linked = {link.target for link in native.links if link.relation == 'supports_rna_sequence'}
    assert {p.key for p in placements} <= linked
    path = tmp_path / 'rna.proteome.json.gz'
    native.save(path)
    assert Proteome.load(path).identity == native.identity


def test_old_sv_export_keeps_unknown_placement_support_unlinked():
    export = export_sv_rna_orfs(reconstruction())
    for occurrence in export['candidates'][0]['occurrences']:
        del occurrence['rna_support']
    proteome = from_rna_result(export, reference=SV_REFERENCE)
    occurrence, = records_of(proteome, 'occurrence')
    assert 'rna_support' not in occurrence.payload
    assert not any(link.target == occurrence.key and link.relation == 'supports_rna_sequence'
                   for link in proteome.links)
    assert len(proteome) == 1


def test_zero_candidate_limited_sv_result_retains_the_event():
    result = reconstruction()
    result['paths'] = []
    result['limitations'] = ['max_paths']
    proteome = from_rna_result(result, reference=SV_REFERENCE)
    assert not len(proteome)
    event = next(r for r in records_of(proteome, 'event') if r.schema_version == 'sv_rna_orfs.v4')
    assert event.status == 'limited' and event.complete is None
    assert event.payload['limitations'] == ('max_paths',)


def test_sv_source_alternatives_and_equal_sequences_do_not_collapse():
    first = reconstruction()
    second = deepcopy(first)
    second['source'] = 'different-library'
    a = from_rna_result(first, reference=SV_REFERENCE)
    b = from_rna_result(second, reference=SV_REFERENCE)
    joint = a.merge(b)
    assert len(joint) == 2
    assert len(joint.by_sequence()['MKMK']) == 2
    assert {p.producer.run_id for p in joint} == {'source', 'different-library'}
    assert a.merge(a).identity == a.identity


def test_fusion_partial_and_unresolved_outcomes_keep_their_provenance(tmp_path):
    fusion, refs, reads = example(partial=1)
    result = reconstruct_fusion(fusion, refs, reads)
    before = deepcopy(result)
    proteome = from_rna_result(result, reference=FUSION_REFERENCE)
    protein, = proteome
    assert protein.sequence == result['paths'][0]['translations'][0]['amino_acids']
    assert protein.producer.start_complete is False
    assert protein.producer.stop_complete is True
    assert protein.producer.origin == 'fusion' and not protein.evidence
    assert set(protein.producer.transcript_ids) == {'donor.1', 'acceptor.1'}
    assert records_of(proteome, 'evidence')[0].to_dict()['payload'] == json.loads(json.dumps(result['evidence']))
    assert result == before
    path = tmp_path / 'fusion.json.gz'
    proteome.save(path)
    assert Proteome.load(path).identity == proteome.identity
    unresolved = reconstruct_fusion(fusion, (), reads)
    empty = from_rna_result(unresolved, reference=FUSION_REFERENCE)
    assert not len(empty)
    assert records_of(empty, 'event')[0].status == 'unresolved'
    assert records_of(empty, 'path')[0].payload['sequence'] == fusion.sequence
    assert records_of(empty, 'evidence')


def test_variant_orf_native_snapshot_keeps_models_support_and_counterfactual(noncoding_rna, tmp_path):
    _, model, observations = noncoding_rna
    variant, evidence = observations(split=True)
    result = assess(variant, evidence, [model])
    before = deepcopy(result)
    proteome = from_rna_result(result, reference=model.context)
    assert len(proteome) == len(result['protein_candidates'])
    assert [p.sequence for p in proteome] == [r['protein_sequence'] for r in result['protein_candidates']]
    assert all(p.producer.origin == 'nc_orf' and not p.evidence for p in proteome)
    rows = records_of(proteome, 'hypothesis')
    assert {r.payload['candidate_id']: r.to_dict()['payload'] for r in rows} == {
        r['candidate_id']: r for r in result['protein_candidates']}
    assert all(not r.payload['initiation_observed'] and not r.payload['translation_observed'] for r in rows)
    assert records_of(proteome, 'model')[0].to_dict()['payload'] == json.loads(json.dumps(result['models'][0]))
    assert from_rna_result(json.loads(json.dumps(result)), reference=model.context).identity == proteome.identity
    assert len(records_of(proteome, 'evidence')) == len(result['evidence_sets'])
    assert result == before
    path = tmp_path / 'nc.proteome.json.gz'
    proteome.save(path)
    assert Proteome.load(path).identity == proteome.identity
    with pytest.raises(ValueError, match='reference context'):
        from_rna_result(result, reference=replace(model.context, annotation_version='another'))


def test_unresolved_ncorf_retains_rna_and_model_outcomes(noncoding_rna):
    _, model, observations = noncoding_rna
    variant, evidence = observations(index=15, alt='C')
    result = assess(variant, evidence, [model])
    proteome = from_rna_result(result, reference=model.context)
    assert not len(proteome)
    assert records_of(proteome, 'model')[0].payload['reason'] == 'reference_start_affected'
    assert records_of(proteome, 'rna_candidate')


def test_dispatch_rejects_unknown_schema_build_and_corrupt_translation():
    assert from_rna_result(exported(), reference=REFERENCE)
    with pytest.raises(ValueError, match='Unsupported'):
        from_rna_result({'schema': 'not-isovar'}, reference=REFERENCE)
    with pytest.raises(ValueError, match='conflicts'):
        from_rna_result(reconstruction(), reference=REFERENCE)
    export = export_sv_rna_orfs(reconstruction())
    export['candidates'][0]['amino_acids'] = 'MMMM'
    with pytest.raises(ValueError, match='translation disagree'):
        from_rna_result(export, reference=SV_REFERENCE)


def test_native_sv_keeps_annotated_frame_products_and_original_paths(tmp_path):
    from tests.test_sv_rna import Scenario, write_bam

    scenario = Scenario()
    result = scenario.run(write_bam(tmp_path / 'rna.bam', scenario.reads()))
    reference = ReferenceContext(9606, 'human', 'test', 'synthetic', '1')
    proteome = from_rna_result(result, reference=reference)
    expected = [t['amino_acids'] for p in result['paths'] for t in p['translations']]
    assert expected
    assert [p.sequence for p in proteome if p.producer.origin == 'sv_translation'] == expected
    assert any(p.producer.origin == 'sv_orf' for p in proteome)
    records = {r.payload['path_id']: r.to_dict()['payload'] for r in records_of(proteome, 'path')}
    for path in result['paths']:
        assert records[path['path_id']] == json.loads(json.dumps({k: v for k, v in path.items() if k != 'translations'}))
    assert any(link.relation == 'has_exploratory_orf_export' for link in proteome.links)
    assert from_rna_result(json.loads(json.dumps(result)), reference=reference).identity == proteome.identity


def test_cli_combines_validated_sources_and_keeps_empty_results(tmp_path):
    reference = tmp_path / 'reference.json'
    reference.write_text(json.dumps(SV_REFERENCE.to_dict()))
    paths, expected = [], None
    for source in ['library-A', 'library-B']:
        raw = reconstruction()
        raw['source'] = source
        value = export_sv_rna_orfs(raw)
        path = tmp_path / (source + '.json.gz')
        path.write_bytes(gzip.compress(json.dumps(value).encode()))
        paths.append(str(path))
        converted = from_rna_result(value, reference=SV_REFERENCE)
        expected = converted if expected is None else expected.merge(converted)
    output, fasta = tmp_path / 'combined.json.gz', tmp_path / 'combined.fasta'
    commands.run(['proteome', '--input', paths[0], '--input', paths[1],
                  '--reference-context', str(reference), '--output', str(output), '--fasta', str(fasta)])
    assert Proteome.load(output).identity == expected.identity
    assert fasta.read_text().count('>') == 2
    empty = export_sv_rna_orfs(dict(reconstruction(), paths=[]))
    empty_path = tmp_path / 'empty.json'
    empty_path.write_text(json.dumps(empty))
    commands.run(['proteome', '--input', str(empty_path), '--reference-context', str(reference),
                  '--output', str(output), '--fasta', str(fasta)])
    restored = Proteome.load(output)
    assert len(restored) == 0 and records_of(restored, 'event')
    assert fasta.read_text() == ''


def test_cli_rejects_collisions_and_reference_mismatch_before_writing(tmp_path, capsys):
    input_path = tmp_path / 'input.json'
    input_path.write_text(json.dumps(export_sv_rna_orfs(reconstruction())))
    reference = tmp_path / 'reference.json'
    reference.write_text(json.dumps(REFERENCE.to_dict()))
    before = input_path.read_bytes()
    with pytest.raises(SystemExit):
        commands.run(['proteome', '--input', str(input_path), '--reference-context', str(reference),
                      '--output', str(input_path)])
    assert 'distinct' in capsys.readouterr().err
    assert input_path.read_bytes() == before
    output = tmp_path / 'output.json.gz'
    with pytest.raises(SystemExit):
        commands.run(['proteome', '--input', str(input_path), '--reference-context', str(reference),
                      '--output', str(output)])
    assert 'conflicts' in capsys.readouterr().err and not output.exists()


def test_fusion_evidence_cannot_be_relabelled_or_inflated():
    result = reconstruct_fusion(*example())
    result['sample_id'] = 'different'
    with pytest.raises(ValueError, match='evidence scope'):
        from_rna_result(result, reference=FUSION_REFERENCE)
    result['provenance']['sample_id'] = 'different'
    with pytest.raises(ValueError, match='evidence scope'):
        from_rna_result(result, reference=FUSION_REFERENCE)
    result = reconstruct_fusion(*example())
    result['evidence']['support']['fragments'] += 1
    with pytest.raises(ValueError, match='support disagrees'):
        from_rna_result(result, reference=FUSION_REFERENCE)


def test_sv_occurrence_support_must_be_part_of_the_union():
    from isovar.rna_evidence import evidence_set

    export = export_sv_rna_orfs(reconstruction())
    candidate = export['candidates'][0]
    candidate['occurrences'][0]['rna_support'].update(evidence_set(['sample', 'source'], [('rg', 'unrelated', 0)]))
    with pytest.raises(ValueError, match='outside the candidate union'):
        from_rna_result(export, reference=SV_REFERENCE)


def test_hcc1395_original_rna_fusion_preserves_complete_sequence_and_support():
    from tests.test_hcc1395_eif3k_fusion import DATA
    from reports.hcc1395_colo829.build_hcc1395_eif3k_fusion import report

    result = report(DATA / 'original_unmapped.bam', DATA / 'models.json.gz',
                    mapped_tiling_bam=DATA / 'coding_tiling.bam')['result']
    reference = ReferenceContext(9606, 'human', 'GRCh37-lite', 'INTEGRATE 2014 hg19 Ensembl exons')
    proteome = from_rna_result(result, reference=reference)
    protein, = proteome
    assert len(protein.sequence) == 531
    assert sha256(protein.sequence.encode()).hexdigest() == (
        '6d2fbb1872c08a507f56a977fb97328471b90a0f077211892d0a5d5ef3053573')
    assert protein.producer.start_complete and protein.producer.stop_complete
    assert not protein.evidence
    path, = records_of(proteome, 'path')
    assert path.payload['junctions'][0]['direct_support']['fragments'] == 149
    assert records_of(proteome, 'translation')[0].payload['translation_observed'] is False


def test_native_sv_unknown_cell_labels_remain_unknown_and_conflicts_fail():
    from isovar.cell_umi import missing_label_row

    result = reconstruction()
    row = missing_label_row(result['observations']['full']['identity'])
    result['cell_umi_evidence']['reads'] = [row]
    proteome = from_rna_result(result, reference=SV_REFERENCE)
    hypothesis, = records_of(proteome, 'hypothesis')
    assert hypothesis.payload['rna_support']['unknown_library_reads'] == 1
    assert hypothesis.payload['rna_support']['umis_complete'] is False
    row['scope'] = {'sample_id': 'other', 'source': 'source'}
    with pytest.raises(ValueError, match='cell/UMI evidence scope'):
        from_rna_result(result, reference=SV_REFERENCE)


def test_candidate_union_cannot_masquerade_as_measured_occurrence_support():
    export = export_sv_rna_orfs(reconstruction())
    candidate = export['candidates'][0]
    candidate['occurrences'][0]['rna_support'] = deepcopy(candidate['rna_support'])
    with pytest.raises(ValueError, match='measured for that occurrence'):
        from_rna_result(export, reference=SV_REFERENCE)


def test_colo829_sv3_keeps_exploratory_alternatives_and_unresolved_frame():
    from reports.hcc1395_colo829.build_unified_sv_manifest import SV3
    from reports.hcc1395_colo829.build_unified_sv_report import snapshots

    snapshot = next(value for _, value in snapshots(event_id=SV3)
                    if value['source_id'] == 'SRR8615617-CNIH3' and value['view_id'] == '1-forward')
    export = export_sv_rna_orfs(snapshot['reconstruction'])
    proteome = from_rna_result(export, reference=REFERENCE)
    assert len(proteome) == len(export['candidates']) > 0
    assert {p.sequence for p in proteome} == {c['amino_acids'] for c in export['candidates']}
    assert all(not p.evidence for p in proteome)
    assert all(r.status == 'exploratory' and not r.payload['translation_observed']
               for r in records_of(proteome, 'hypothesis'))
    occurrences = records_of(proteome, 'occurrence')
    target = [r for r in occurrences if r.payload['path_id'] == '653d5d9c015d3e590fb149b3']
    assert target and all(r.payload['frame_status'] == 'no_gene_anchor' for r in target)
    assert records_of(proteome, 'event')[0].status == 'limited'
