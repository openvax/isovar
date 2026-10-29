"""Source preparation must retain naming, orientation and annotation uncertainty."""

from types import SimpleNamespace

from pyensembl import EnsemblRelease
import pytest
from varcode import ContigResolver
from varcode.structural_variant import StructuralVariant

from isovar import prepare_sv_source_record, read_sv_source
from isovar._file_identity import file_sha256


def source(tmp_path, contig='virus_A', mate='decoy_B', alt=None):
    path = tmp_path / 'source.vcf'
    path.write_text(f'{contig}\t100\tevent\tN\t{alt or "N[" + mate + ":200["}\t.\tPASS\tCIPOS=-2,3\n')
    return read_sv_source(path, sha256=file_sha256(path), reference_name='GRCh38')[0]


def test_missing_annotation_retains_custom_contig_rna_views(tmp_path):
    record = source(tmp_path)
    prepared = prepare_sv_source_record(record, EnsemblRelease(87),
        ContigResolver(['1'], reference_name='GRCh38'), annotation='Ensembl 87', region_flank=0)
    assert prepared.references == prepared.hypotheses == ()
    assert len(prepared.limitations) == 3
    forward, reverse = prepared.views
    assert forward.donor.contig == 'virus_A'
    assert (forward.donor.position, forward.donor.strand) == (100, '+')
    assert (forward.acceptor.position, forward.acceptor.strand) == (199, '+')
    assert (reverse.donor.position, reverse.donor.strand) == (199, '-')
    assert (reverse.acceptor.position, reverse.acceptor.strand) == (100, '-')
    assert forward.regions == (('virus_A', 99, 100), ('decoy_B', 199, 200))
    assert prepared.nominal_igv == ('virus_A:100-100', 'decoy_B:200-200')
    assert record.variant.ci_start == (-2, 3)


def test_annotation_aliases_do_not_rewrite_source_models(tmp_path, monkeypatch):
    record = source(tmp_path, contig='chr1', mate='chr1')
    genome = EnsemblRelease(87)
    transcript = SimpleNamespace(id='tx1', sequence='ATGAAATAA', strand='+',
        exons=[SimpleNamespace(start=97, end=105)], contains_start_codon=True, contains_stop_codon=True,
        start_codon_complete=True, stop_codon_complete=True,
        first_start_codon_spliced_offset=0, last_stop_codon_spliced_offset=8)
    gene = SimpleNamespace(name='GENE1', transcripts=[transcript])
    queried = []

    def genes_at_locus(contig, start, end):
        queried.append(contig)
        return [gene]

    monkeypatch.setattr(genome, 'genes_at_locus', genes_at_locus)
    monkeypatch.setattr(StructuralVariant, 'effects', lambda self: [])
    prepared = prepare_sv_source_record(record, genome, ContigResolver(['1'], reference_name='GRCh38'),
                                         annotation='test annotation')
    assert set(queried) == {'1'}
    assert len(prepared.references) == 1
    reference, = prepared.references
    assert reference.contig == 'chr1' and reference.cds_start == 0 and reference.cds_end == 9
    assert prepared.reference_payload()['models'][0]['contig'] == 'chr1'
    assert prepared.gene_names == ('GENE1',)
    assert not prepared.limitations


@pytest.mark.parametrize('field,value', [('region_flank', -1), ('annotation_flank', True), ('annotation', '')])
def test_invalid_preparation_context(tmp_path, field, value):
    kwargs = dict(annotation='test')
    kwargs[field] = value
    with pytest.raises(ValueError):
        prepare_sv_source_record(source(tmp_path), EnsemblRelease(87),
                                 ContigResolver(['1'], reference_name='GRCh38'), **kwargs)
