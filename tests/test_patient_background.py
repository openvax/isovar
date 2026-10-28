"""Patient baselines use exact cDNA edits and separate original REF evidence."""

from copy import deepcopy
from dataclasses import replace
import gzip
import json

import pytest

from isovar import AlleleRead, IsovarResult, ProteinSequence, export_patient_background_pairs
from isovar.cli import commands
from isovar.dna import reverse_complement_dna
from isovar.genetic_code import translate_cdna
from isovar.known_variants import KnownVariants
from isovar.read_evidence import ReadEvidence
from isovar.variant_sequence import VariantSequence
from tests.test_edit_attribution import (
    CHANGING_GERMLINE, CO_SOMATIC, FOCAL, SILENT_GERMLINE, TRANSCRIPT, translation, variant,
)
from tests.test_edit_attribution import ntf3 as ntf3


def read(prefix, allele, suffix, name='read', strand='+', transcript_ids=None):
    if strand == '-':
        prefix, allele, suffix = (reverse_complement_dna(suffix), reverse_complement_dna(allele),
                                  reverse_complement_dna(prefix))
    return AlleleRead(prefix, allele, suffix, name,
                      compatible_transcript_ids=transcript_ids,
                      source_alignments=[(('rg', name, 0), None)])


def result(t=None, *, known=None, ref_reads=(), reads=None):
    if t is None:
        t = translation('GCATGGCC', 'AAGATTAA', 'MAQD')
    if known is None:
        known = KnownVariants(somatic=[t.reference_context.variant, CO_SOMATIC],
                              germline=[SILENT_GERMLINE, CHANGING_GERMLINE])
    seq = t.untrimmed_variant_sequence
    if reads is None:
        reads = [read(seq.prefix, seq.alt, seq.suffix, strand=t.reference_context.strand)]
    t.untrimmed_variant_sequence = VariantSequence(seq.prefix, seq.alt, seq.suffix, reads)
    return IsovarResult(
        variant=t.reference_context.variant, predicted_effect=None,
        sorted_protein_sequences=[ProteinSequence.from_translations([t], known)],
        read_evidence=ReadEvidence(108, t.reference_context.variant.ref, t.reference_context.variant.alt,
                                   ref_reads=list(ref_reads), alt_reads=reads, other_reads=[]))


def export(value, **kwargs):
    return export_patient_background_pairs([value], sample_id='patient', source='RNA', **kwargs)


def pair(value, **kwargs):
    row, = export(value, **kwargs)['events'][0]['pairs']
    return row


def test_focal_pair_keeps_germline_cosomatic_and_unexplained_bases():
    value = result()
    p = pair(value)
    assert p['state'] == 'paired'
    assert p['mutant']['protein_sequence'] == 'MAQD'
    baseline = p['baseline']
    assert baseline['protein_sequence'] == 'MAKD'  # annotation is MAKG
    assert baseline['nucleotide_sequence'] == 'GCATGGCCAAAGATTAA'
    assert baseline['construction'] == 'transcript_cdna_somatic_reversion'
    assert baseline['observation'] == 'computed_only'
    assert not baseline['normal_tissue_observed'] and not p['translation_observed']
    assert p['phase']['state'] == 'single_original_read_span'
    assert [e['origin'] for e in p['removed_edits']] == ['nominated_variant']
    assert {e['origin'] for e in p['retained_edits']} == {'germline_variant', 'co_somatic_variant', 'unexplained'}
    both = pair(value, baseline='supplied_somatic_haplotype')
    assert both['baseline']['nucleotide_sequence'] == 'GCATGGCCAAAGACTAA'
    assert {e['origin'] for e in both['removed_edits']} == {'nominated_variant', 'co_somatic_variant'}
    assert all(e['origin'] != 'germline_variant' for e in both['removed_edits'])


def test_reference_allele_is_observed_only_with_separate_exact_spans():
    baseline_read = read('GCATGGCC', 'A', 'AAGATTAA', name='ref')
    value = result(ref_reads=[baseline_read, read('CCATGGCT', 'A', 'AAGGCTAA', name='annotation')])
    p = pair(value)
    b = p['baseline']
    assert b['observation'] == 'reference_allele_rna_span'
    support = b['reference_allele_coding_support']['single_read_span_support']
    assert support['fragments'] == 1
    payload = export(value)
    assert payload['evidence_sets'][support['evidence_set_id']]['evidence_scope'] == ['patient', 'RNA']
    assert not b['normal_tissue_observed']
    assert pair(value, baseline='supplied_somatic_haplotype')['baseline']['observation'] == 'computed_only'
    wrong_path = read('GCATGGCC', 'A', 'AAGATTAA', transcript_ids={'different'})
    assert pair(result(ref_reads=[wrong_path]))['baseline']['observation'] == 'computed_only'


@pytest.mark.parametrize('ref,alt,suffix', [('A', '', 'AAGACTAA'), ('AAA', '', 'GACTAA'), ('AAA', 'CC', 'GACTAA')])
def test_focal_frameshift_inframe_and_complex_indels_retranslate_cdna(ref, alt, suffix):
    focal = variant(108, ref, alt)
    reference_suffix = 'CCATGGCTAAAGGCTAA'[8 + len(ref):]
    t = translation('CCATGGCC', suffix, '', focal=focal, reference_suffix=reference_suffix,
                    context_after=reference_suffix)
    t.amino_acids, t.ends_with_stop_codon = translate_cdna(t.cdna_sequence[2:])
    value = result(t, known=KnownVariants(somatic=[focal], germline=[SILENT_GERMLINE, CHANGING_GERMLINE]))
    p = pair(value)
    assert p['state'] == 'paired'
    assert p['baseline']['nucleotide_sequence'] == 'CCATGGCCAAAGACTAA'
    assert p['baseline']['protein_sequence'] == 'MAKD'
    if ref == 'A':
        assert not p['baseline']['ends_with_stop_codon']
        assert p['baseline']['trailing_incomplete_codon_bases'] == 1
    if ref == 'AAA' and alt == '':
        assert p['mutant']['protein_sequence'] == 'MAD'
        assert p['baseline']['ends_with_stop_codon']


def test_germline_deletion_remains_in_counterfactual():
    t = translation('CCATGGCT', 'GAAAC', '', reference_suffix='GAAAAC', context_after='GAAAAC')
    t.amino_acids, t.ends_with_stop_codon = translate_cdna(t.cdna_sequence[2:])
    p = pair(result(t, known=KnownVariants(somatic=[FOCAL], germline=[variant(110, 'A', '')])))
    assert p['state'] == 'paired'
    assert p['baseline']['nucleotide_sequence'] == 'CCATGGCTAGAAAC'
    assert [e['origin'] for e in p['retained_edits']] == ['germline_variant']
    assert p['retained_edits'][0]['alt_bases'] == ''


def test_minus_strand_reference_reads_and_edit_coordinates():
    class MinusTranscript:
        id, name, gene = 'minus', 'Minus', TRANSCRIPT.gene
        contig, start, end, strand = '1', 100, 116, '-'

        def spliced_offset(self, position):
            return self.end - position

    t = translation('CCATGGCC', 'AAGACTAA', 'MAQD')
    context = t.reference_context
    context.strand, context.variant, context.transcripts = '-', variant(108, 'T', 'G'), (MinusTranscript(),)
    known = KnownVariants(somatic=[context.variant], germline=[variant(109, 'A', 'G'), variant(104, 'C', 'T')])
    p = pair(result(t, known=known, ref_reads=[read('CCATGGCC', 'A', 'AAGACTAA', 'ref', strand='-')]))
    assert p['state'] == 'paired' and p['baseline']['protein_sequence'] == 'MAKD'
    assert p['baseline']['observation'] == 'reference_allele_rna_span'
    assert p['removed_edits'][0]['cdna_interval'] == [8, 9]
    assert p['removed_edits'][0]['reference_bases'] == 'A'
    assert {e['origin'] for e in p['retained_edits']} == {'germline_variant'}


def test_original_sid_ntf3_retains_supplied_proximal_germline_in_baseline(ntf3):
    focal, adjacent, run = ntf3
    original, = run(variants=[focal], germline_variants=[adjacent])
    payload = export(original)
    rows = [r for r in payload['events'][0]['pairs'] if r['state'] == 'paired']
    assert rows and all(r['removed_edits'][0]['origin'] == 'nominated_variant' for r in rows)
    assert all(any(e['origin'] == 'germline_variant' for e in r['retained_edits']) for r in rows)
    assert all(r['mutant']['rna_support']['fragments'] > 0 for r in rows)
    assert all(r['baseline']['construction'] == 'transcript_cdna_somatic_reversion' for r in rows)
    assert any(r['baseline']['protein_sequence'] != r['mutant']['protein_sequence'] for r in rows)


def test_tiled_and_missing_phase_are_explicit_and_partial_ref_is_not_a_span():
    reads = [read('GCATGGCC', 'C', 'AAG', 'left'), read('GCC', 'C', 'AAGATTAA', 'right')]
    p = pair(result(reads=reads, ref_reads=[read('GCATGGCC', 'A', 'AAG', 'ref')]))
    assert p['phase']['state'] == 'assembly_only'
    assert p['baseline']['observation'] == 'computed_only'
    assert p['baseline']['reference_allele_coding_support']['single_read_span_support']['reads'] == 0
    assert pair(result(reads=[]))['phase']['state'] == 'incomplete_rna_coverage'


def test_alternative_transcript_offsets_and_incomplete_replay_are_not_collapsed():
    t = translation('GCATGGCC', 'AAGATTAA', 'MAQD')
    other = replace(TRANSCRIPT, id='tx-other', start=99)
    t.reference_context.transcripts = (TRANSCRIPT, other)
    rows = export(result(t))['events'][0]['pairs']
    assert len(rows) == 2 and len({r['pair_id'] for r in rows}) == 2
    assert {r['reference']['cdna_interval'][0] for r in rows} == {0, 1}
    assert {r['baseline']['protein_sequence'] for r in rows} == {'MAKD'}
    t = translation('GCATGGCC', 'AAGATTAA', 'MAQD', reference_suffix='AAGGCTAAGGG')
    bounded = pair(result(t))
    assert bounded['state'] == 'paired'
    assert bounded['reference']['nucleotide_sequence'] == 'CCATGGCTAAAGGCTAA'
    assert bounded['baseline']['nucleotide_sequence'] == 'GCATGGCCAAAGATTAA'
    t = translation('GCATGGCC', 'AAGGCTAAGGG', 'MAQG')
    p = pair(result(t))
    assert p['state'] == 'unavailable' and p['baseline'] is None
    assert p['reason'] == 'transcript_edits_do_not_reproduce_observed_rna'


def test_capped_frame_boundary_has_no_fabricated_baseline():
    focal = variant(108, 'AAA', 'CCCC')
    t = translation('CCATGGCT', 'GGCTAA', '', focal=focal, reference_suffix='GGCTAA', context_after='GGCTAA')
    amino, _ = translate_cdna(t.cdna_sequence[2:])
    t.amino_acids, t.ends_with_stop_codon = amino[:3], False
    p = pair(result(t, known=KnownVariants(somatic=[focal])))
    assert p['state'] == 'unavailable' and p['reason'] == 'focal_allele_outside_paired_coding_window'


def test_empty_and_scoped_export_is_stable_without_mutating_input():
    value = result()
    before = deepcopy(value)
    first = export(value)
    assert first == export(value) and value == before
    other = export_patient_background_pairs([value], sample_id='other', source='RNA')
    def nonempty(x):
        return {k for k, row in x['evidence_sets'].items() if row['reads']}
    assert nonempty(first).isdisjoint(nonempty(other))
    empty = value.clone_with_updates(sorted_protein_sequences=[])
    event, = export(empty)['events']
    assert not event['pairs'] and event['state'] == 'no_protein_hypotheses'
    with pytest.raises(ValueError, match='baseline'):
        export(value, baseline='normal')


def test_cli_uses_same_export_and_rejects_input_overwrite(tmp_path, monkeypatch, capsys):
    from isovar.cli import isovar_patient_background_pairs as cli
    value = result()
    monkeypatch.setattr(cli, 'run_isovar_from_parsed_args', lambda options: [value])
    output = tmp_path / 'pairs.json.gz'
    args = ['patient-background-pairs', '--vcf', 'somatic.vcf', '--bam', 'rna.bam',
            '--sample-id', 'patient', '--source', 'RNA', '--output', str(output)]
    commands.run(args)
    with gzip.open(output, 'rt') as handle:
        assert json.load(handle) == export(value)
    original = output.read_bytes()
    with pytest.raises(SystemExit):
        commands.run(args + ['--vcf', str(output)])
    assert output.read_bytes() == original and 'overwrite' in capsys.readouterr().err
