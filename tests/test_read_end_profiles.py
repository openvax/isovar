"""Keep sequence sources, processing facts and cached matches reproducible."""

from dataclasses import asdict, replace
import json

import pytest

from isovar import (Adapter, ReadEndProfile, ReadLibraryProvenance, ReadProcessingStep,
                    ReadEndInferenceCache, infer_read_ends, read_end_profile_from_fasta)
from isovar._file_identity import file_sha256
from isovar.read_end_profiles import read_end_profiles_from_json


MOTIF = 'AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT'
BODY = 'GCTACGTCGCTG'
PROFILE = ReadEndProfile('explicit', adapters=(Adapter('motif', MOTIF),))


def test_explicit_provenance_round_trip_and_fingerprint():
    provenance = ReadLibraryProvenance(platform='PACBIO', library_state='segmented',
        processing=(ReadProcessingStep('skera', 'split', version='1.4.0', source='run-receipt.json'),),
        sources=('sample-sheet.tsv',))
    profile = replace(PROFILE, provenance=provenance)
    assert profile.provenance.kit is None and profile.provenance.instrument is None
    assert ReadEndProfile.from_dict(asdict(profile)) == profile
    assert profile.fingerprint != PROFILE.fingerprint
    assert replace(profile, provenance=replace(provenance, kit='supplied-kit')).fingerprint != profile.fingerprint


@pytest.mark.parametrize('role', ['adapter', 'primer', 'linker'])
def test_roles_preserved_in_both_orientations(role):
    from isovar.dna import reverse_complement_dna
    profile = replace(PROFILE, adapters=(replace(PROFILE.adapters[0], role=role, source='sequence-source'),))
    for reverse in (False, True):
        sequence = reverse_complement_dna(BODY + MOTIF) if reverse else BODY + MOTIF
        view = infer_read_ends(sequence, profile=profile, reverse=reverse, trim_adapters=True)
        assert view.sequence == (reverse_complement_dna(BODY) if reverse else BODY)
        assert {annotation.kind for annotation in view.annotations} == {role}


def test_fasta_api_and_relative_json_share_the_same_profile(tmp_path):
    fasta = tmp_path / 'primers.fa'
    fasta.write_text('>first comment\n' + MOTIF + '\n>second\nGCTAGCATGCTACGTACGATCGTA\n')
    digest = file_sha256(fasta)
    fields = dict(name='future-library', version='1', both_orientations=True,
                  provenance=dict(preparation='custom', sources=['sample-sheet.tsv']))
    defaults = dict(role='primer', end='left', mate=2, min_overlap=12, max_error_rate=0.1)
    expected = read_end_profile_from_fasta(fasta, sha256=digest, profile=fields, adapter_defaults=defaults)
    path = tmp_path / 'profile.json'
    path.write_text(json.dumps(dict(fields, adapter_fasta=dict(path='primers.fa', sha256=digest, defaults=defaults))))
    assert read_end_profiles_from_json(path) == expected
    assert [adapter.name for adapter in expected.adapters] == ['first', 'second']
    assert all(adapter.source == 'sha256:' + digest for adapter in expected.adapters)
    path.write_text(json.dumps(dict(read_groups={'group1': json.loads(path.read_text())})))
    assert read_end_profiles_from_json(path) == {'group1': expected}
    with pytest.raises(ValueError, match='checksum mismatch'):
        read_end_profile_from_fasta(fasta, sha256='0' * 64, profile=fields)


@pytest.mark.parametrize('contents', ['', '>same\nACGTACGTACGT\n>same\nACGTACGTACGT\n',
    '@first\nACGTACGTACGT\n+\nFFFFFFFFFFFF\n', '>ambiguous\nACGTNNNNACGT\n'])
def test_fasta_rejects_empty_duplicate_fastq_and_ambiguous_sequences(tmp_path, contents):
    path = tmp_path / 'adapters.fa'
    path.write_text(contents)
    with pytest.raises(ValueError):
        read_end_profile_from_fasta(path, sha256=file_sha256(path), profile=dict(name='custom'))


def test_cache_reuses_matches_but_not_quality_placement_or_tail_trimming():
    sequence = BODY + 'A' * 25 + MOTIF
    cache = ReadEndInferenceCache()
    first = infer_read_ends(sequence, [30] * len(sequence), PROFILE, trim_adapters=True,
                            trim_poly_a=True, cache=cache)
    second = infer_read_ends(sequence, None, PROFILE, soft_clip_bounds=(0, len(sequence)),
                             trim_adapters=True, trim_poly_a=True, upstream_poly_a_length=-1, cache=cache)
    assert first.sequence == BODY
    assert second.sequence == sequence and second.qualities is None
    assert second.upstream_poly_a_length == -1
    assert cache.hits == 1 and cache.misses == 1
    assert first.original_qualities == (30,) * len(sequence)
    assert first.annotations == second.annotations
    assert second == infer_read_ends(sequence, None, PROFILE, soft_clip_bounds=(0, len(sequence)),
                                     trim_adapters=True, trim_poly_a=True, upstream_poly_a_length=-1)


def test_cache_separates_profile_mate_and_orientation_and_bounds_memory():
    sequence = BODY + MOTIF
    cache = ReadEndInferenceCache(max_entries=2, max_bases=2 * len(sequence))
    for kwargs in ({}, dict(mate=2), dict(reverse=True), dict(profile=replace(PROFILE, version='2'))):
        infer_read_ends(sequence, cache=cache, **dict(dict(profile=PROFILE), **kwargs))
    assert cache.misses == 4 and cache.hits == 0
    assert cache.retained_bases <= 2 * len(sequence)
    infer_read_ends(sequence, profile=PROFILE, cache=cache)
    assert cache.misses == 5
    oversized = sequence * 3
    infer_read_ends(oversized, cache=cache)
    assert cache.retained_bases <= 2 * len(sequence)


@pytest.mark.parametrize('kwargs', [dict(mate=True), dict(mate=3), dict(reverse='yes'),
    dict(hard_clips=(-1, 0)), dict(hard_clips=(0,)), dict(hard_clips=(0, False))])
def test_inference_rejects_invalid_orientation_and_unavailable_lengths(kwargs):
    with pytest.raises(ValueError):
        infer_read_ends(BODY, **kwargs)


def test_collector_settings_accept_typed_inline_and_pinned_profiles(tmp_path):
    from isovar import ReadCollector
    definition = asdict(PROFILE)
    inline = ReadCollector.from_dict(dict(read_end_profile=definition, trim_adapters=True))
    typed = ReadCollector.from_dict(dict(read_end_profile=PROFILE, trim_adapters=True))
    path = tmp_path / 'profile.json'
    path.write_text(json.dumps(definition))
    pinned = ReadCollector.from_dict(dict(read_end_profile={'file_id': 'profile'}, trim_adapters=True),
                                    files={'profile': path})
    assert inline.read_end_profile == typed.read_end_profile == pinned.read_end_profile == PROFILE
    assert definition == asdict(PROFILE)
    grouped = ReadCollector.from_dict(dict(read_end_profile={'read_groups': {'rg': definition}}))
    assert grouped.read_end_profile == {'rg': PROFILE}
    with pytest.raises(ValueError, match='undeclared file'):
        ReadCollector.from_dict(dict(read_end_profile={'file_id': 'missing'}))
    with pytest.raises(ValueError, match='cannot be combined'):
        ReadCollector.from_dict(dict(read_end_profile={'file_id': 'profile', 'adapters': []}))


def test_sv_manifest_source_profile_and_external_fasta_are_checked_on_resume(tmp_path, monkeypatch):
    from isovar import sv_rna_audit
    from tests.test_sv_rna_audit import manifest_fixture
    manifest_path, manifest = manifest_fixture(tmp_path)
    fasta = tmp_path / 'adapters.fa'
    fasta.write_text('>adapter\n' + MOTIF + '\n')
    profile_path = tmp_path / 'profile.json'
    profile_path.write_text(json.dumps(dict(name='source-specific',
        adapter_fasta=dict(path=fasta.name, sha256=file_sha256(fasta)))))
    manifest['files']['profile'] = dict(path=profile_path.name, sha256=file_sha256(profile_path))
    manifest['read_collection'] = dict(read_end_profile=dict(name='default-library'))
    manifest['sources'][0]['read_collection'] = dict(read_end_profile={'file_id': 'profile'})
    sv_rna_audit.write_json(manifest_path, manifest)
    seen = []
    original = sv_rna_audit.prepare_rna_library

    def prepare(*args, **kwargs):
        seen.append(kwargs['read_collector'].read_end_profile.name)
        return original(*args, **kwargs)

    monkeypatch.setattr(sv_rna_audit, 'prepare_rna_library', prepare)
    output = tmp_path / 'output'
    first = sv_rna_audit.run_sv_rna_audit(manifest_path, output)
    assert seen == ['source-specific']
    assert first == sv_rna_audit.run_sv_rna_audit(manifest_path, output, resume=True)
    fasta.write_text('>adapter\n' + MOTIF[:-1] + 'A\n')
    with pytest.raises(ValueError, match='FASTA checksum mismatch'):
        sv_rna_audit.run_sv_rna_audit(manifest_path, output, resume=True)


def test_cohort_profile_file_is_part_of_acquisition_and_scientific_identity(tmp_path):
    from isovar.cohort import _selected_file_keys, _small_variant_file_keys, run_small_variant_source
    from isovar.sv_rna_audit import write_json
    from tests.testing_helpers import data_path
    from pathlib import Path
    calls = tmp_path / 'calls.json'
    write_json(calls, [dict(contig='14', start=92811641, ref='GAGA', alt='G')])
    profile = tmp_path / 'profile.json'
    profile.write_text(json.dumps(dict(name='explicit-unknown', adapters=[])))
    bam = Path(data_path('data/colo829_whole_run/SRR8615617.all_sites.grch38.bam'))
    task = dict(source_id='rna', variants='calls', ensembl_release=87,
                read_collection=dict(read_end_profile={'file_id': 'profile'}))
    source = dict(source_id='rna', sample_id='COLO829', bam='bam')
    files = dict(calls=calls, bam=bam, profile=profile)
    manifest = dict(files={key: {} for key in files}, sources=[source], small_variants=[task])
    assert set(_small_variant_file_keys(task, source)) == set(files)
    assert _selected_file_keys(manifest, None) == set(files)
    receipt = run_small_variant_source(task, source, files, tmp_path / 'cohort')
    assert receipt['status'] == 'completed'
    assert receipt['identity']['inputs']['profile'] == file_sha256(profile)
