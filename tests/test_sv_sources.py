"""Source interpretation must preserve defects and never manufacture support."""

from collections import Counter
from pathlib import Path

import pytest

from isovar.cli.commands import run
from isovar.sv_rna_audit import file_sha256, read_json
from isovar.sv_sources import read_sv_source, sv_source_report


DATA = Path(__file__).parent / 'data'
COLO = DATA / 'colo829_whole_run/colo829-sv-grch38.vcf'
HCC = DATA / 'sv_sources/HCC1395.svscope.pcr.vcf'


def inspect(path, profile='standard', **kwargs):
    return read_sv_source(path, sha256=file_sha256(path), reference_name='GRCh38', profile=profile, **kwargs)


def test_colo_profile_resolves_only_documented_breakend_encoding():
    strict = {record.record_id: record for record in inspect(COLO)}
    records = inspect(COLO, 'colo829_paulin_paz_v2')
    assert len(records) == 54
    assert Counter(r.state for r in records) == dict(resolved_adjacency=46, invalid=1, unresolved=7)
    repaired = [r for r in records if r.transformations]
    assert {r.record_id for r in repaired} == {
        'COLO829_SOMATIC_SV22_BND', 'COLO829_SOMATIC_SV28_BND',
        'COLO829_SOMATIC_SV29_BND', 'COLO829_SOMATIC_SV45_BND'}
    for record in repaired:
        assert record.original_row == strict[record.record_id].original_row
        assert strict[record.record_id].state == 'invalid'
        assert record.state == 'resolved_adjacency'
    first = repaired[0].variant.junctions[0]
    assert [(e.contig, e.position, e.keeps) for e in first] == [
        ('chr7', 151049571, 'left'), ('chr15', 84141972, 'right')]
    inversion = next(r for r in records if r.record_id.endswith('SV25_INV'))
    assert inversion.variant is None and 'END greater than POS' in inversion.reasons[0]
    # A benchmark nomination's supplied 0/0 genotype must not be rewritten.
    assert next(r for r in repaired if 'SV29' in r.record_id).original_row.endswith('\tGT\t0/0')


def test_hcc_pcr_flags_and_unresolved_coordinate_claims_are_separate():
    records = inspect(HCC, 'hcc1395_svscope_2026')
    assert len(records) == 1118
    assert sum(r.pcr_validation is True for r in records) == 32
    assert sum(r.pcr_validation is None for r in records) == 1086
    assert Counter(r.state for r in records) == dict(resolved_adjacency=157, invalid=780, unresolved=181)
    inversion = next(r for r in records if r.record_id == 'Minda_192')
    assert inversion.variant is None
    assert {c.mate_position for c in inversion.coordinate_claims} == {169157265}
    assert all(c.mate_position - c.position == 11174717 for c in inversion.coordinate_claims)
    assert 'SVLEN=11174710' in inversion.original_row
    assert len([r for r in records if r.reasons == ('source_consensus_alleles_require_reference_anchoring',)]) == 22
    assert all(r.variant is None for r in records if r.record_id.startswith('TDscope.'))


def test_named_source_profile_rejects_modified_source_bytes(tmp_path):
    path = tmp_path / 'changed.vcf'
    path.write_bytes(COLO.read_bytes() + b'\n')
    with pytest.raises(ValueError, match='pinned GRCh38 original'):
        inspect(path, 'colo829_paulin_paz_v2')
    with pytest.raises(ValueError, match='checksum mismatch'):
        read_sv_source(path, sha256='0' * 64, reference_name='GRCh38')


@pytest.mark.parametrize('info,expected', [('PCR_Validation=True', True),
    ('PCR_Validation=False', False), ('.', None)])
def test_pcr_validation_has_three_states(tmp_path, info, expected):
    path = tmp_path / 'PCR_Validated.vcf'
    path.write_text('virusA\t100\tevent\tN\tN[decoy_B:200[\t.\tPASS\t' + info + '\n')
    record, = inspect(path)
    assert record.pcr_validation is expected
    assert record.state == 'resolved_adjacency'
    assert [end.contig for end in record.variant.junctions[0]] == ['virusA', 'decoy_B']


@pytest.mark.parametrize('info', ['END=100;END=200', 'PCR_Validation=maybe', 'CIPOS=not,integers',
    'END=200;SVTYPE', 'END=200.5', 'END=200;CIPOS=1,-1', 'END=200;CIPOS=-1,.'])
def test_invalid_fields_are_retained_without_a_variant(tmp_path, info):
    path = tmp_path / 'bad.vcf'
    line = 'chr1\t100\tevent\tN\t<DEL>\t.\tPASS\t' + info
    path.write_text(line + '\n')
    record, = inspect(path)
    assert record.state == 'invalid' and record.variant is None
    assert record.original_row == line


def test_source_cli_matches_api(tmp_path):
    path = tmp_path / 'report.json.gz'
    run(['sv-source', '--vcf', str(HCC), '--sha256', file_sha256(HCC),
         '--reference-name', 'GRCh38', '--profile', 'hcc1395_svscope_2026', '--output', str(path)])
    expected = sv_source_report(inspect(HCC, 'hcc1395_svscope_2026'), sha256=file_sha256(HCC),
                                reference_name='GRCh38', profile='hcc1395_svscope_2026')
    assert read_json(path) == expected


def test_annotation_adapter_resolves_names_without_mutating_source():
    from pyensembl import EnsemblRelease
    from varcode import ContigResolver
    record = next(r for r in inspect(COLO, 'colo829_paulin_paz_v2') if r.record_id.endswith('SV22_BND'))
    genome = EnsemblRelease(87)
    adapted = record.for_annotation(genome, ContigResolver(['7', '15'], reference_name='GRCh38'))
    assert (adapted.contig, adapted.mate_contig) == ('7', '15')
    assert (record.variant.contig, record.variant.mate_contig) == ('chr7', 'chr15')
    with pytest.raises(ValueError, match='Unresolved annotation contig'):
        record.for_annotation(genome, ContigResolver(['7'], reference_name='GRCh38'))
    with pytest.raises(ValueError, match='assembly mismatch'):
        record.for_annotation(genome, ContigResolver(['7', '15'], reference_name='GRCh37'))


def test_breakend_confidence_intervals_survive_shared_parser(tmp_path):
    path = tmp_path / 'breakend.vcf'
    path.write_text('chr1\t100\tevent\tN\tN[chr2:200[\t.\tPASS\tCIPOS=-5,7;CIEND=-8,9\n')
    record, = inspect(path)
    assert record.state == 'resolved_adjacency'
    assert record.variant.ci_start == (-5, 7)
    assert record.variant.ci_end == (-8, 9)
