"""One resolved code follows native RNA candidates into every consumer (#700)."""

from copy import deepcopy
import json

import pytest
from orfdata import TranscriptORF
from pyensembl import Genome
from pyensembl.species import find_species_by_name
from varcode import Variant

from isovar import AlleleRead, ReadEvidence, from_rna_result
from isovar.dna import reverse_complement_dna
from .genomes_for_testing import grch38
from .test_variant_rna_orfs import assess
from .test_variant_rna_orfs import noncoding_rna as noncoding_rna
from .test_explicit_genetic_codes import native_orf as native_orf


def observations(reference, offset=201):
    transcript = reference.transcript
    positions = [p for start, end in sorted(transcript.exon_intervals)
                 for p in range(start, end + 1)]
    reverse = transcript.strand == '-'
    if reverse:
        positions.reverse()
    index = reference.start + offset
    sequence = reference.transcript_sequence
    ref, alt = sequence[index], 'C' if sequence[index] != 'C' else 'A'
    prefix, suffix = sequence[max(0, index - 72):index], sequence[index + 1:index + 94]
    if reverse:
        ref, alt = reverse_complement_dna(ref), reverse_complement_dna(alt)
        prefix, suffix = reverse_complement_dna(suffix), reverse_complement_dna(prefix)
    variant = Variant(transcript.contig, positions[index], ref, alt, transcript.genome)
    evidence = ReadEvidence.from_variant_and_allele_reads(
        variant, [AlleleRead(prefix, alt, suffix, name='synthetic-code-control')])
    return variant, evidence


@pytest.mark.parametrize('name,table', [('MT-ND5-201', 2), ('MT-ND6-201', 2), ('TP53-001', 1)])
def test_annotated_code_matches_explicit_control_and_native_replay(name, table, tmp_path):
    from orfdata import Proteome

    native = TranscriptORF.from_annotation(grch38.transcripts_by_name(name)[0])
    assert native.translation_table is None
    explicit = TranscriptORF.at(
        native.transcript, native.start, end=native.end, orf_id='explicit-control',
        source='synthetic-control', protein=native.protein_sequence, translation_table=table)
    variant, evidence = observations(native)
    result = assess(variant, evidence, [native])
    control = assess(variant, evidence, [explicit])
    assert result['protein_candidates']
    fields = ('protein_sequence', 'reference_protein_window',
              'reference_allele_counterfactual_protein', 'ends_with_stop_codon')
    assert {tuple(row[k] for k in fields) for row in result['protein_candidates']} == {
        tuple(row[k] for k in fields) for row in control['protein_candidates']}
    for candidate in result['protein_candidates']:
        assert candidate['translation_table'] == table
        assert candidate['translation_table_source'] == 'reference_context'
        assert candidate['rna_support']['fragments'] == 1
        start = candidate['reference_protein_start_0based']
        window = candidate['reference_protein_window']
        assert window == native.protein_sequence[start:start + len(window)]
        assert candidate['reference_allele_counterfactual_protein'] == window
    assert all(row['translation_table_source'] == 'orf_model' for row in control['protein_candidates'])
    assert result['models'][0]['reference']['translation_table'] is None
    assert result['allele_support']['fragments'] == 1
    assert result['schema'] == 'isovar.variant_rna_orfs.v2'
    replay = json.loads(json.dumps(result))
    proteome = from_rna_result(replay, reference=native.context)
    assert {p.sequence for p in proteome} == {c['protein_sequence'] for c in result['protein_candidates']}
    path = tmp_path / 'native-code.proteome.json.gz'
    proteome.save(path)
    assert Proteome.load(path).identity == proteome.identity


@pytest.mark.parametrize('strand', ['+', '-'])
def test_nuclear_numt_like_annotation_uses_nuclear_code(tmp_path, strand):
    # Deliberately synthetic nuclear carrier: mitochondrial-like sequence and
    # name are not evidence of mitochondrial location or translation.
    coding = 'ATG' + 'GCC' * 4 + 'ATA' + 'GCC' * 4 + 'TGA' + 'GCC' * 4 + 'TAA'
    sequence = 'C' * 15 + coding + 'G' * 15
    attrs = ('gene_id "NUC"; gene_name "MT-ND5"; transcript_id "NUC-T"; '
             'transcript_name "MT-ND5-NUMT"; gene_biotype "protein_coding"; '
             'transcript_biotype "protein_coding"; protein_id "NUC-P";')
    rows = []
    for kind, start, end in [('gene', 0, len(sequence)), ('transcript', 0, len(sequence)),
                             ('exon', 0, len(sequence)), ('CDS', 15, 45),
                             ('start_codon', 15, 18), ('stop_codon', 45, 48)]:
        low, high = ((100 + start, 99 + end) if strand == '+' else
                     (100 + len(sequence) - end, 99 + len(sequence) - start))
        extra = ' exon_id "E"; exon_number "1";' if kind == 'exon' else ''
        rows.append(f'1\tsynthetic\t{kind}\t{low}\t{high}\t.\t{strand}\t0\t{attrs}{extra}\n')
    gtf, cdna, protein = [tmp_path / name for name in ('reference.gtf', 'cdna.fa', 'protein.fa')]
    gtf.write_text(''.join(rows))
    cdna.write_text('>NUC-T\n' + sequence + '\n')
    protein.write_text('>NUC-P\nMAAAAIAAAA\n')
    genome = Genome('synthetic-numt', 'synthetic', '1', gtf_path_or_url=str(gtf),
                    transcript_fasta_paths_or_urls=[str(cdna)], protein_fasta_paths_or_urls=[str(protein)],
                    cache_directory_path=str(tmp_path))
    genome.species = find_species_by_name('homo_sapiens')
    genome.index()
    try:
        native = TranscriptORF.from_annotation(genome.transcript_by_id('NUC-T'))
        assert native.translation_table is None and native.protein_sequence == 'MAAAAIAAAA'
        result = assess(*observations(native, 18), [native])
        assert result['protein_candidates']
        for candidate in result['protein_candidates']:
            assert candidate['translation_table'] == 1
            assert candidate['translation_table_source'] == 'reference_context'
            assert candidate['reference_protein_window'] == 'MAAAAIAAAA'
            assert candidate['ends_with_stop_codon']
    finally:
        genome.close()


def test_v1_replay_retains_undeclared_source_code(noncoding_rna):
    _, model, observe = noncoding_rna
    variant, evidence = observe()
    current = assess(variant, evidence, [model])
    historical = deepcopy(current)
    historical['schema'] = 'isovar.variant_rna_orfs.v1'
    historical['models'][0]['reference']['translation_table'] = None
    for row in historical['protein_candidates']:
        row.pop('translation_table', None)
        row.pop('translation_table_source', None)
    snapshot = from_rna_result(historical, reference=model.context)
    assert len(snapshot) == len(current['protein_candidates'])
    assert historical['models'][0]['reference']['translation_table'] is None


def test_native_source_defined_codes_survive_reconstruction(native_orf):
    _, reference = native_orf
    offset = min(201, 3 * (len(reference.protein_sequence) // 2))
    result = assess(*observations(reference, offset), [reference])
    assert result['protein_candidates']
    for candidate in result['protein_candidates']:
        assert candidate['translation_table'] == reference.translation_table
        assert candidate['translation_table_source'] == 'orf_model'
        start = candidate['reference_protein_start_0based']
        window = candidate['reference_protein_window']
        assert window == reference.protein_sequence[start:start + len(window)]
    assert len(from_rna_result(result, reference=reference.context)) == len(result['protein_candidates'])


@pytest.mark.parametrize('field,value', [
    ('translation_table', None), ('translation_table', True), ('translation_table', '1'),
    ('translation_table', 999), ('translation_table', 27), ('translation_table', 2),
    ('translation_table_source', None), ('translation_table_source', 'reference_context'),
])
def test_v2_snapshot_rejects_invalid_or_conflicting_effective_code(noncoding_rna, field, value):
    _, model, observe = noncoding_rna
    variant, evidence = observe()
    result = assess(variant, evidence, [model])
    result['protein_candidates'][0][field] = value
    with pytest.raises(ValueError, match='translation.table|recoding|provenance'):
        from_rna_result(result, reference=model.context)


def test_v2_table_reader_never_defaults_missing_fields():
    from isovar import variant_orf_translation_table

    source = {'translation_table': None}
    with pytest.raises(ValueError, match='requires an integer'):
        variant_orf_translation_table({}, source, schema='isovar.variant_rna_orfs.v2')
    assert variant_orf_translation_table({}, source, schema='isovar.variant_rna_orfs.v1') is None
    with pytest.raises(ValueError, match='Unsupported native'):
        variant_orf_translation_table({}, source, schema='isovar.variant_rna_orfs.v0')
