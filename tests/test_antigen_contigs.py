"""Antigen references use the same exact/alias/ambiguity contract as variants."""

from hashlib import md5
from types import SimpleNamespace

import pysam
import pytest
from varcode import ContigResolver

from isovar.antigen_reference import contig_in_reference, validate_alignment_reference
from isovar.contig_resolution import UnresolvedContigError


@pytest.mark.parametrize('query,names,expected', [
    ('MT', ['CHRM'], 'CHRM'), ('chrM', ['MT'], 'MT'),
    ('chr1', ['1'], '1'), ('chr1', ['1', 'chr1'], 'chr1'),
    ('HPV16', ['HPV16', 'HPV18'], 'HPV16'), ('chrEBV', ['chrEBV'], 'chrEBV'),
    ('hs38d1', ['hs38d1'], 'hs38d1'), ('0', ['0', 'HPV16'], '0'),
    ('GL000207.1', ['GL000207.1', 'GL000207.2'], 'GL000207.1'),
    ('chrUn_KI270442v1', ['chrUn_KI270442v1'], 'chrUn_KI270442v1'),
])
def test_antigen_reference_exact_and_unique_aliases(query, names, expected):
    assert contig_in_reference(query, names) == expected
    assert contig_in_reference(query, ContigResolver(names)) == expected


@pytest.mark.parametrize('query,names,status', [
    ('chr1', ['1', 'CHR1'], 'ambiguous'),
    ('HPV16', ['HPV18'], 'missing'), ('GL000207.1', ['GL000207.2'], 'missing'),
])
def test_unresolved_names_keep_structured_diagnostics(query, names, status):
    with pytest.raises(UnresolvedContigError) as caught:
        contig_in_reference(query, names)
    assert caught.value.resolution.status == status
    assert caught.value.to_dict()['original_name'] == query
    if status == 'ambiguous':
        assert set(caught.value.resolution.candidates) == set(names)


def alignment_reference(fasta, rows):
    """Construct in-memory sequence/header inputs without annotation fixtures."""
    genome = SimpleNamespace(fasta=fasta, reference_name='GRCh38',
                             annotation_name='test', annotation_version='1')
    bam = SimpleNamespace(header=pysam.AlignmentHeader.from_dict({'SQ': rows}))
    return bam, genome


@pytest.mark.parametrize('bam_name,fasta_name', [
    ('1', 'chr1'), ('chr1', '1'), ('chrM', 'MT'),
    ('KI270728.1', 'chr16_KI270728v1_random'),
    ('chrUn_KI270442v1', 'KI270442.1'),
    ('HPV16', 'HPV16'), ('chrEBV', 'chrEBV'), ('hs38d1', 'hs38d1'),
])
def test_alignment_reference_aliases_retain_names_and_validate_sequence(bam_name, fasta_name):
    sequence = 'ACGT'
    bam, genome = alignment_reference({fasta_name: sequence}, [
        dict(SN=bam_name, LN=4, M5=md5(sequence.encode()).hexdigest())])
    result = validate_alignment_reference(bam, genome, contigs=[bam_name])
    assert result['checked_sequence_contigs'] == [bam_name]
    assert result['bam_sequence_identity_unverified'] == []
    assert result['alignment_contig_resolutions'] == [dict(
        bam_contig=bam_name, fasta_contig=fasta_name,
        status='exact' if bam_name == fasta_name else 'alias')]


def test_alignment_reference_exact_name_wins_with_different_alias_sequence():
    bam, genome = alignment_reference({'chr1': 'AAAA', '1': 'TTTT'}, [
        dict(SN='1', LN=4, M5=md5(b'TTTT').hexdigest())])
    assert validate_alignment_reference(bam, genome, contigs=['1'])['checked_sequence_contigs'] == ['1']


def test_alignment_reference_alias_without_checksum_remains_unverified():
    bam, genome = alignment_reference({'MT': 'ACGT'}, [dict(SN='chrM', LN=4)])
    result = validate_alignment_reference(bam, genome, contigs=['MT', 'chrM'])
    assert result['checked_sequence_contigs'] == []
    assert result['bam_sequence_identity_unverified'] == ['chrM']


@pytest.mark.parametrize('fasta,status', [
    ({'1': 'ACGT', 'CHR1': 'ACGT'}, 'ambiguous'),
    ({'2': 'ACGT'}, 'missing'),
])
def test_alignment_reference_unresolved_names_keep_diagnostics(fasta, status):
    bam, genome = alignment_reference(fasta, [dict(SN='chr1', LN=4)])
    with pytest.raises(UnresolvedContigError) as caught:
        validate_alignment_reference(bam, genome)
    assert caught.value.resolution.status == status


@pytest.mark.parametrize('row,match', [
    (dict(SN='KI270728.1', LN=5), 'length mismatch'),
    (dict(SN='KI270728.1', LN=4, M5='0' * 32), 'checksum mismatch'),
])
def test_alignment_reference_alias_does_not_hide_sequence_mismatch(row, match):
    bam, genome = alignment_reference({'chr16_KI270728v1_random': 'ACGT'}, [row])
    with pytest.raises(ValueError, match=match):
        validate_alignment_reference(bam, genome, contigs=['KI270728.1'])


def test_distinct_bam_contigs_cannot_collapse_onto_one_fasta_sequence():
    bam, genome = alignment_reference({'MT': 'ACGT'}, [
        dict(SN='MT', LN=4), dict(SN='chrM', LN=4)])
    with pytest.raises(ValueError, match='Distinct BAM contigs'):
        validate_alignment_reference(bam, genome)
