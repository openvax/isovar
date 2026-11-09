"""An unusable insertion/model pair must not abort the other supplied models."""

import pytest
from orfdata import ORFModel, ORFRecord
from pyensembl import Genome
from pyensembl.species import find_species_by_name
from varcode import Variant

from isovar import AlleleRead, ProteinSequenceCreator, ReadEvidence, reconstruct_variant_orfs
from isovar.dna import reverse_complement_dna
from isovar.reference_sequence_key import ReferenceSequenceKey
from isovar.variant_helpers import NoExonicReferenceBases, interbase_range_affected_by_variant_on_transcript


@pytest.fixture(params=['+', '-'])
def insertion_models(tmp_path, request):
    strand = request.param
    sequence = 'C' * 15 + 'ATG' + 'GCC' * 10 + 'TAA' + 'G' * 15
    exons = {'valid': [(100, 165)],
             'intronic': [(100, 129), (145, 165)] if strand == '+' else [(136, 165), (100, 120)]}
    sequences = {'valid': sequence, 'intronic': sequence[:30] + sequence[45:]}
    rows = []
    for name, intervals in exons.items():
        attrs = (f'gene_id "G"; gene_name "GENE"; transcript_id "{name}"; '
                 'gene_biotype "lncRNA"; transcript_biotype "lncRNA";')
        rows.extend(f'1\ttest\t{kind}\t100\t165\t.\t{strand}\t.\t{attrs}\n'
                    for kind in (('gene', 'transcript') if name == 'valid' else ('transcript',)))
        rows.extend(f'1\ttest\texon\t{start}\t{end}\t.\t{strand}\t.\t{attrs} '
                    f'exon_id "{name}-{i}"; exon_number "{i + 1}";\n'
                    for i, (start, end) in enumerate(intervals))
    gtf, fasta = tmp_path / 'reference.gtf', tmp_path / 'reference.fa'
    gtf.write_text(''.join(rows))
    fasta.write_text(''.join(f'>{name}\n{seq}\n' for name, seq in sequences.items()))
    genome = Genome('insertion_test', 'synthetic', '1', gtf_path_or_url=str(gtf),
                    transcript_fasta_paths_or_urls=[str(fasta)], cache_directory_path=str(tmp_path))
    genome.species = find_species_by_name('homo_sapiens')
    genome.index()
    references = []
    for name, intervals in exons.items():
        blocks = tuple(sorted((max(start, 115), min(end, 150)) for start, end in intervals))
        model = ORFModel(name, name, 'G', 'GENE', '1', strand, blocks,
                         'M' + 'A' * (10 if name == 'valid' else 5), 'lncRNA',
                         'insertion_test', 9606, 'homo_sapiens', source='synthetic')
        references.append(ORFRecord.from_model(model).bind(genome.transcript_by_id(name)))
    genomic = sequence if strand == '+' else reverse_complement_dna(sequence)
    inserted = 'GCC' if strand == '+' else 'GGC'
    variant = Variant('1', 132, '', inserted, genome)
    reads = [AlleleRead(genomic[:33], inserted, genomic[33:], name=f'read-{i}',
                        reference_blocks=[(0, 33, 99, 132), (36, 69, 132, 165)]) for i in range(2)]
    try:
        yield genome, references, variant, reads
    finally:
        genome.close()


@pytest.mark.parametrize('anchor', [99, 132, 166])
def test_nonexonic_insertion_anchor_has_typed_outcome(insertion_models, anchor):
    genome, (_, intronic), _, _ = insertion_models
    variant = Variant('1', anchor, '', 'GCC', genome)
    with pytest.raises(NoExonicReferenceBases, match='insertion anchor') as caught:
        interbase_range_affected_by_variant_on_transcript(variant, intronic.transcript)
    assert type(caught.value.__cause__) is ValueError
    assert ReferenceSequenceKey.from_variant_and_transcript(variant, intronic.transcript, 30) is None


def test_exonic_insertion_offsets_keep_strand_and_boundary_conventions(insertion_models):
    genome, (_, intronic), _, _ = insertion_models
    expected = ([(100, 1), (129, 30), (145, 31), (165, 51)] if intronic.strand == '+' else
                [(165, 0), (136, 29), (120, 30), (100, 50)])
    for anchor, offset in expected:
        variant = Variant('1', anchor, '', 'GCC', genome)
        assert interbase_range_affected_by_variant_on_transcript(variant, intronic.transcript) == (offset, offset)


@pytest.mark.parametrize('reverse_models', [False, True])
@pytest.mark.parametrize('with_reads', [False, True])
def test_intronic_model_remains_unresolved_without_aborting_valid_model(insertion_models, reverse_models, with_reads):
    _, references, variant, reads = insertion_models
    evidence = ReadEvidence.from_variant_and_allele_reads(variant, reads if with_reads else [])
    result = reconstruct_variant_orfs(
        variant, evidence, references[::-1] if reverse_models else references,
        sample_id='sample', source='original-library',
        protein_sequence_creator=ProteinSequenceCreator(min_variant_sequence_coverage=1))
    models = {row['reference']['resolved_transcript']['id']: row for row in result['models']}
    assert models['intronic']['state'] == 'unresolved'
    assert models['intronic']['reason'] == 'no_exonic_reference_bases'
    assert models['intronic']['candidate_ids'] == []
    assert result['allele_support']['reads'] == (2 if with_reads else 0)
    if with_reads:
        assert models['valid']['state'] == 'rna_frame_candidates'
        assert all(candidate['transcript_id'] == 'valid' for candidate in result['protein_candidates'])
        candidate = next(candidate for candidate in result['protein_candidates']
                         if candidate['protein_sequence'] == 'M' + 'A' * 11)
        assert candidate['rna_support']['fragments'] == 2
    else:
        assert models['valid']['reason'] == 'no_alt_reads'
        assert result['protein_candidates'] == []


def test_unexpected_anchor_lookup_errors_still_propagate(insertion_models, monkeypatch):
    _, (_, intronic), variant, _ = insertion_models

    def broken_lookup(position):
        raise RuntimeError('broken annotation storage')

    monkeypatch.setattr(intronic.transcript, 'spliced_offset', broken_lookup)
    with pytest.raises(RuntimeError, match='broken annotation storage'):
        interbase_range_affected_by_variant_on_transcript(variant, intronic.transcript)
