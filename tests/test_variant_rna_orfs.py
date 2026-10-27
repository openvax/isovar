"""Alternate coding frames retain actual RNA and never inherit translation evidence."""

from dataclasses import replace
import json
from pathlib import Path

import pytest
from orfdata import ORFModel, ORFRecord
from pyensembl import Genome
from pyensembl.species import find_species_by_name
from varcode import Variant

from isovar import AlleleRead, ProteinSequenceCreator, ReadEvidence, reconstruct_variant_orfs
from isovar.dna import reverse_complement_dna
from isovar.reference_context_helpers import reference_contexts_for_variant


@pytest.fixture(params=['+', '-'])
def noncoding_rna(tmp_path, request):
    strand = request.param
    sequence = 'C' * 15 + 'ATG' + 'GCC' * 10 + 'TAA' + 'G' * 15
    end = 99 + len(sequence)
    attrs = 'gene_id "G"; gene_name "NONCODING"; transcript_id "T"; gene_biotype "lncRNA"; transcript_biotype "lncRNA";'
    gtf, fasta = tmp_path / 'reference.gtf', tmp_path / 'reference.fa'
    gtf.write_text(''.join(f'1\ttest\t{kind}\t100\t{end}\t.\t{strand}\t.\t{attrs}{extra}\n'
                           for kind, extra in [('gene', ''), ('transcript', ''),
                                                ('exon', ' exon_id "E"; exon_number "1";')]))
    fasta.write_text('>T\n' + sequence + '\n')
    genome = Genome('test_assembly', 'test_annotation', '1', gtf_path_or_url=str(gtf),
                    transcript_fasta_paths_or_urls=[str(fasta)], cache_directory_path=str(tmp_path))
    genome.species = find_species_by_name('homo_sapiens')
    genome.index()
    transcript = genome.transcript_by_id('T')
    start, stop = (115, 150) if strand == '+' else (end - 50, end - 15)
    record = ORFRecord.from_model(ORFModel(
        'noncoding_orf', 'T', 'G', 'NONCODING', '1', strand, ((start, stop),),
        'M' + 'A' * 10, 'lncRNA', 'test_assembly', 9606, 'homo_sapiens', source='synthetic'))
    reference = record.bind(transcript)

    def observations(index=34, alt='T', split=False):
        genomic = sequence if strand == '+' else reverse_complement_dna(sequence)
        position = index if strand == '+' else len(sequence) - index - 1
        allele = alt if strand == '+' else reverse_complement_dna(alt)
        variant = Variant('1', 100 + position, genomic[position], allele, genome)
        spans = [(0, 50), (20, len(sequence))] if split else [(0, len(sequence))] * 2
        reads = []
        for i, (left, right) in enumerate(spans):
            assert left <= position < right
            reads.append(AlleleRead(
                genomic[left:position], allele, genomic[position + 1:right], name=f'fragment{i}',
                reference_blocks=[(0, right - left, 99 + left, 99 + right)]))
        return variant, ReadEvidence.from_variant_and_allele_reads(variant, reads)

    yield genome, reference, observations
    genome.close()


def assess(variant, evidence, references):
    return reconstruct_variant_orfs(
        variant, evidence, references, sample_id='sample', source='RNA-library',
        protein_sequence_creator=ProteinSequenceCreator(
            protein_sequence_length=49, protein_sequence_preference='support',
            min_variant_sequence_coverage=1))


def test_noncoding_parent_recovers_observed_mutant_protein(noncoding_rna):
    _, reference, observations = noncoding_rna
    variant, evidence = observations()
    assert reference_contexts_for_variant(variant) == []
    result = assess(variant, evidence, [reference])
    candidate, = result['protein_candidates']
    assert candidate['protein_sequence'] == 'MAAAAAVAAAA'
    assert candidate['reference_protein_window'] == 'MAAAAAAAAAA'
    assert candidate['reference_allele_counterfactual_protein'] == 'MAAAAAAAAAA'
    assert candidate['mutation_interval'] == [6, 7]
    assert candidate['focal_amino_acid_change']
    assert candidate['sequence_extent'] == 'complete_coding_candidate'
    assert candidate['coding_support']['single_read_span_support']['fragments'] == 2
    assert candidate['coding_support']['minimum_fragment_depth'] == 2
    assert not candidate['initiation_observed'] and not candidate['translation_observed']
    assert result['models'][0]['reference']['resolved_transcript']['biotype'] == 'lncRNA'
    assert candidate['transcript_id'] == 'T'
    assert json.loads(json.dumps(result))['schema'] == 'isovar.variant_rna_orfs.v1'


def test_synonymous_context_is_not_promoted_to_amino_acid_change(noncoding_rna):
    _, reference, observations = noncoding_rna
    variant, evidence = observations(index=35)
    candidate, = assess(variant, evidence, [reference])['protein_candidates']
    assert candidate['protein_sequence'] == 'MAAAAAAAAAA'
    assert not candidate['focal_amino_acid_change']


def test_tiling_does_not_count_as_single_read_span(noncoding_rna):
    _, reference, observations = noncoding_rna
    variant, evidence = observations(split=True)
    result = assess(variant, evidence, [reference])
    full = [p for p in result['protein_candidates'] if p['sequence_extent'] == 'complete_coding_candidate']
    assert full
    for candidate in full:
        assert candidate['coding_support']['fully_tiled']
        assert candidate['coding_support']['single_read_span_support']['fragments'] == 0
        assert min(candidate['coding_support']['fragment_depth']) == 1


def test_no_model_keeps_rna_but_does_not_guess_a_frame(noncoding_rna):
    _, _, observations = noncoding_rna
    variant, evidence = observations()
    result = assess(variant, evidence, [])
    assert result['assembly_state'] == 'assembled'
    assert result['rna_candidates']
    assert not result['models'] and not result['protein_candidates']
    assert result['reference_incompatible_support']['fragments'] == 2


def test_no_alt_and_start_loss_are_explicit_unresolved_outcomes(noncoding_rna):
    _, reference, observations = noncoding_rna
    variant, evidence = observations()
    empty = ReadEvidence.from_variant_and_allele_reads(variant, [])
    result = assess(variant, empty, [reference])
    assert result['assembly_state'] == 'no_alt_reads'
    assert result['models'][0]['reason'] == 'no_alt_reads'
    variant, evidence = observations(index=15, alt='C')
    result = assess(variant, evidence, [reference])
    assert result['models'][0]['reason'] == 'reference_start_affected'
    assert not result['protein_candidates']


def test_reference_and_focal_evidence_validation(noncoding_rna):
    genome, reference, observations = noncoding_rna
    variant, evidence = observations()
    with pytest.raises(ValueError, match='Duplicate supplied ORF'):
        assess(variant, evidence, [reference, reference])
    different = Variant('1', variant.start + 1, variant.ref, variant.alt, genome)
    with pytest.raises(ValueError, match='Read evidence does not describe'):
        assess(different, evidence, [reference])
    mismatched = replace(reference.context, annotation_version='other')
    old = reference.context
    object.__setattr__(reference, 'context', mismatched)
    try:
        with pytest.raises(ValueError, match='reference contexts differ'):
            assess(variant, evidence, [reference])
    finally:
        object.__setattr__(reference, 'context', old)


@pytest.mark.parametrize('removed,inserted,frameshift', [(3, '', False), (1, '', True), (0, 'GCC', False)])
def test_indels_keep_observed_flanks_and_coding_support(noncoding_rna, removed, inserted, frameshift):
    genome, reference, _ = noncoding_rna
    seq, index = reference.transcript_sequence, 33
    genomic = seq if reference.strand == '+' else reverse_complement_dna(seq)
    pos = index if reference.strand == '+' else len(seq) - index - removed
    allele = inserted if reference.strand == '+' else reverse_complement_dna(inserted)
    # Unanchored literal alleles use Varcode's interbase insertion convention.
    variant = Variant('1', 100 + pos - int(not removed), genomic[pos:pos + removed], allele, genome)
    reads = [AlleleRead(genomic[:pos], allele, genomic[pos + removed:], name=f'indel{i}') for i in range(2)]
    evidence = ReadEvidence.from_variant_and_allele_reads(variant, reads)
    candidate, = assess(variant, evidence, [reference])['protein_candidates']
    assert candidate['frameshift'] == frameshift
    assert candidate['focal_amino_acid_change']
    assert candidate['coding_support']['minimum_fragment_depth'] == 2
    assert candidate['coding_support']['single_read_span_support']['fragments'] == 2
    assert candidate['reference_allele_counterfactual_protein'] == reference.protein_sequence
    if not frameshift:
        assert len(candidate['protein_sequence']) == len(reference.protein_sequence) + (len(inserted) - removed) // 3


def test_competing_placements_never_become_independent_support(noncoding_rna):
    _, reference, observations = noncoding_rna
    variant, evidence = observations()
    read = evidence.alt_reads[0]
    alternatives = [AlleleRead(read.prefix, read.allele, read.suffix, name='same',
                              source_alignments=((('', 'same', 0), (0, start, '66M', False)),))
                    for start in [99, 100]]
    evidence = ReadEvidence.from_variant_and_allele_reads(variant, alternatives)
    result = assess(variant, evidence, [reference])
    assert result['allele_support']['fragments'] == 1
    assert len(result['protein_candidates']) == 2
    assert len({p['candidate_id'] for p in result['protein_candidates']}) == 2
    assert all(p['rna_support']['fragments'] == 1 for p in result['protein_candidates'])
    assert all(p['coding_support']['minimum_fragment_depth'] == 1 for p in result['protein_candidates'])


def test_observed_premature_stop_is_retained_with_zero_width_change(noncoding_rna):
    genome, reference, _ = noncoding_rna
    seq, index = reference.transcript_sequence, 33
    genomic = seq if reference.strand == '+' else reverse_complement_dna(seq)
    pos = index if reference.strand == '+' else len(seq) - index - 3
    allele = 'TAA' if reference.strand == '+' else 'TTA'
    variant = Variant('1', 100 + pos, genomic[pos:pos + 3], allele, genome)
    reads = [AlleleRead(genomic[:pos], allele, genomic[pos + 3:], name=f'stop{i}') for i in range(2)]
    result = assess(variant, ReadEvidence.from_variant_and_allele_reads(variant, reads), [reference])
    candidate, = result['protein_candidates']
    assert candidate['protein_sequence'] == 'MAAAAA'
    assert candidate['mutation_interval'] == [6, 6]
    assert candidate['ends_with_stop_codon']
    assert candidate['focal_amino_acid_change']
    assert candidate['coding_interval'][1] - candidate['coding_interval'][0] == 21
    assert candidate['coding_support']['minimum_fragment_depth'] == 2


def test_colo829_dynll2_original_rna_separates_merged_and_single_read_spans():
    import orfdata
    import pysam
    from pyensembl import EnsemblRelease
    from isovar import ReadCollector

    genome = EnsemblRelease(87)
    variant = Variant('17', 58089866, 'G', 'A', genome)
    record = orfdata.load_catalog().get('c17norep133')
    reference = record.bind(genome.transcript_by_id(record.transcript_id))
    path = Path(__file__).parent / 'data/colo829_whole_run/SRR8615617.all_sites.grch38.bam'
    with pysam.AlignmentFile(path) as bam:
        evidence = ReadCollector().read_evidence_for_variant(variant, bam)
    result = reconstruct_variant_orfs(variant, evidence, [reference], sample_id='COLO829', source='SRR8615617')
    assert result['allele_support']['fragments'] == 45
    candidates = [p for p in result['protein_candidates'] if p['protein_sequence'] == 'MPSLEALRCVCLWCV']
    assert candidates
    candidate = max(candidates, key=lambda p: p['rna_support']['fragments'])
    assert candidate['reference_protein_window'] == 'MPSLGALRCVCLWCV'
    assert candidate['mutation_interval'] == [4, 5]
    assert candidate['rna_support']['fragments'] == 30
    assert candidate['coding_support']['allele_observation_span_support']['fragments'] == 21
    assert candidate['coding_support']['single_read_span_support']['fragments'] == 15
    assert candidate['coding_support']['minimum_fragment_depth'] == 24
    assert candidate['sequence_extent'] == 'partial_coding_candidate'
    assert not candidate['translation_observed']
