"""Paired local coding windows with explicit somatic-reversion provenance."""

from . import __version__
from .coding_interval_support import coding_interval_support
from .genetic_code import translate_cdna
from .protein_hypotheses import _EvidenceSets, _check_translation, _completeness, _edits, _event_id, _hypothesis_id, _variant
from .rna_evidence import content_identifier
from .transcript_edit_helpers import transcript_edit_from_variant


SCHEMA = 'isovar.patient_background_pairs.v1'
BASELINES = ('focal', 'supplied_somatic_haplotype')


class _Unavailable(ValueError):
    """A specific RNA/context cannot justify the requested baseline."""


def _replay(translation, transcript, edits):
    orf, context = translation.variant_orf, translation.reference_context
    focal = transcript_edit_from_variant(context.variant, transcript)
    reference = (orf.reference_cdna_sequence_before_variant + context.sequence_at_variant_locus
                 + orf.reference_cdna_sequence_after_variant)
    left = focal.cdna_start - len(orf.reference_cdna_sequence_before_variant)
    right = left + len(reference)
    cursor, observed, positioned = left, '', []
    for edit in edits:
        begin, end = edit['cdna_interval']
        if not left <= begin <= end <= right or begin < cursor:
            raise _Unavailable('overlapping_or_out_of_bounds_transcript_edits')
        observed += reference[cursor - left:begin - left]
        positioned.append(dict(edit, observed_interval=[len(observed), len(observed) + len(edit['alt_bases'])],
                               reference_bases=reference[begin - left:end - left]))
        observed += edit['alt_bases']
        cursor = end
    observed += reference[cursor - left:]
    # Annotation matching commonly supplies a longer downstream flank than
    # the RNA reaches. A terminal reference overhang is not a deletion: trim
    # it only after every retained edit and an exact replay of observed bases.
    extra = len(observed) - len(orf.cdna_sequence)
    if (extra < 0 or observed[:len(orf.cdna_sequence)] != orf.cdna_sequence
            or any(e['observed_interval'][1] > len(orf.cdna_sequence) for e in positioned)):
        raise _Unavailable('transcript_edits_do_not_reproduce_observed_rna')
    if extra:
        reference = reference[:-extra]
        right -= extra
    return positioned, dict(transcript_id=transcript.id, cdna_interval=[left, right],
                            nucleotide_sequence=reference, sequence_id=content_identifier('nt', reference))


def _boundary(position, edits, *, right=False):
    shift = 0
    for edit in edits:
        begin, end = edit['observed_interval']
        if begin < position < end:
            raise _Unavailable('reverted_edit_crosses_coding_boundary')
        if end < position or (end == position and (begin < end or right)):
            shift += len(edit['reference_bases']) - (end - begin)
    return position + shift


def _phase(support):
    if support['single_read_span_support']['reads']:
        return 'single_original_read_span'
    if support['allele_observation_span_support']['reads']:
        return 'allele_observation_span'
    return 'assembly_only' if support['fully_tiled'] else 'incomplete_rna_coverage'


def _pair(translation, transcript, known, result, evidence, baseline):
    _check_translation(translation)
    orf, context = translation.variant_orf, translation.reference_context
    sequence = orf.cdna_sequence
    start = orf.offset_to_first_complete_codon
    end = start + 3 * (len(translation.amino_acids) + int(translation.ends_with_stop_codon))
    def support(reads):
        return evidence.support(reads, 'exact_paired_coding_interval_spans')
    mutant_support = coding_interval_support(
        translation.reads, sequence, orf.variant_cdna_interval_start, (start, end), context.strand, support)
    row = dict(
        transcript_id=transcript.id, strand=context.strand, translation_table=context.translation_table,
        state='unavailable', reason=None, baseline_policy=baseline,
        mutant=dict(nucleotide_sequence=sequence, nucleotide_sequence_id=content_identifier('nt', sequence),
                    coding_interval=[start, end],
                    protein_sequence=translation.amino_acids, ends_with_stop_codon=translation.ends_with_stop_codon,
                    rna_support=evidence.support(translation.reads, 'reads_assembled_into_this_cdna'),
                    coding_support=mutant_support),
        phase=dict(state=_phase(mutant_support), scope='the_mutant_coding_window',
                   normal_tissue_observed=False),
        baseline=None, reference=None, removed_edits=[], retained_edits=[],
        translation_observed=False)
    try:
        selected = sorted((e for e in _edits(translation, known) if e['transcript_id'] == transcript.id),
                          key=lambda e: e['cdna_interval'])
        positioned, row['reference'] = _replay(translation, transcript, selected)
        origins = {'nominated_variant'} if baseline == 'focal' else {'nominated_variant', 'co_somatic_variant'}
        removed = [e for e in positioned if e['origin'] in origins]
        row['removed_edits'] = removed
        row['retained_edits'] = [e for e in positioned if e['origin'] not in origins]
        if sum(e['origin'] == 'nominated_variant' for e in removed) != 1:
            raise _Unavailable('focal_edit_not_uniquely_attributed')
        if not start <= orf.variant_cdna_interval_start <= orf.variant_cdna_interval_end <= end:
            raise _Unavailable('focal_allele_outside_paired_coding_window')
        for edit in removed:
            begin, finish = edit['observed_interval']
            if finish <= start and (len(edit['reference_bases']) - (finish - begin)) % 3:
                raise _Unavailable('removed_frameshift_precedes_coding_window')
        counter = sequence
        for edit in reversed(removed):
            begin, finish = edit['observed_interval']
            counter = counter[:begin] + edit['reference_bases'] + counter[finish:]
        counter_start = _boundary(start, removed)
        counter_end = _boundary(end, removed, right=True)
        dna = counter[counter_start:counter_end]
        amino, stop = translate_cdna(dna, translation_table=context.translation_table)
        translated_end = counter_start + 3 * (len(amino) + int(stop))
        ref_reads = [read for read in result.read_evidence.ref_reads
                     if read.compatible_transcript_ids is None or transcript.id in read.compatible_transcript_ids]
        ref_support = coding_interval_support(
            ref_reads, counter, _boundary(orf.variant_cdna_interval_start, removed),
            (counter_start, translated_end), context.strand, support)
        observed = _phase(ref_support)
        row.update(state='paired' if amino else 'zero_length_baseline',
                   baseline=dict(
                       construction='transcript_cdna_somatic_reversion',
                       observation=('reference_allele_rna_span' if observed in
                                    ('single_original_read_span', 'allele_observation_span') else 'computed_only'),
                       nucleotide_sequence=counter, nucleotide_sequence_id=content_identifier('nt', counter),
                       coding_interval=[counter_start, translated_end],
                       available_coding_interval=[counter_start, counter_end],
                       protein_sequence=amino, ends_with_stop_codon=stop,
                       trailing_incomplete_codon_bases=0 if stop else len(dna) % 3,
                       reference_allele_coding_support=ref_support,
                       reference_allele_phase=observed,
                       normal_tissue_observed=False, translation_observed=False))
    except _Unavailable as error:
        row['reason'] = str(error)
    return row


def export_patient_background_pairs(isovar_results, *, sample_id, source, baseline='focal'):
    """Pair each nucleotide/transcript hypothesis with a declared RNA background.

    Parameters
    ----------
    isovar_results : iterable of IsovarResult
        Native results, preferably from ``run_isovar`` with supplied germline
        calls. All retained proteins and compatible transcript contexts are used.
    sample_id, source : str
        Sample and read-set identities used to scope original RNA evidence.
    baseline : {'focal', 'supplied_somatic_haplotype'}
        Revert only the nominated focal edit (default), or every edit attributed
        to the supplied somatic calls on this RNA. Germline and unexplained
        changes remain. The export lists actual removed and retained edits.

    Returns
    -------
    dict
        JSON-ready ``isovar.patient_background_pairs.v1`` export. Baselines are
        computed, bounded counterfactuals; separate reference-allele RNA spans
        are measured when present. No reference padding or normal-tissue/protein
        expression is inferred. Ambiguous replay/frame contexts remain explicit.
    """
    if baseline not in BASELINES:
        raise ValueError(f'Unknown patient-background baseline: {baseline!r}')
    if any(not isinstance(value, str) or not value.strip() for value in (sample_id, source)):
        raise ValueError('sample_id and source must be nonempty strings')
    evidence = _EvidenceSets([sample_id, source])
    events, supplied_variants = [], {}
    for result in isovar_results:
        event_id = _event_id(result.variant)
        pairs = []
        for protein in result.sorted_protein_sequences:
            hypothesis_id = _hypothesis_id(protein, event_id)
            known = protein.known_variants
            if known is not None:
                for origin in ('somatic', 'germline'):
                    for variant in getattr(known, origin):
                        supplied_variants[_event_id(variant)] = dict(origin=origin, variant=_variant(variant))
            for translation in protein.translations:
                for transcript in sorted(translation.reference_context.transcripts, key=lambda t: t.id):
                    pair = _pair(translation, transcript, protein.known_variants, result, evidence, baseline)
                    pair['hypothesis_id'] = hypothesis_id
                    pair['translation_id'] = content_identifier('translation', [
                        event_id, translation.cdna_sequence, translation.offset_to_first_complete_codon,
                        sorted(t.id for t in translation.reference_context.transcripts),
                        translation.amino_acids, translation.reference_context.translation_table])
                    pair['pair_id'] = content_identifier('patient_background_pair', [event_id, pair])
                    pairs.append(pair)
        pairs.sort(key=lambda row: row['pair_id'])
        limit, complete = _completeness(result)
        events.append(dict(event_id=event_id, variant=_variant(result.variant), pairs=pairs,
                           state='retained' if pairs else 'no_protein_hypotheses',
                           protein_hypotheses_complete=complete, protein_sequence_limit=limit,
                           phased_variants_in_supporting_reads=sorted(
                               _event_id(v) for v in result.phased_variants_in_supporting_reads)))
    return dict(schema=SCHEMA, isovar_version=__version__, sample_id=sample_id, source=source,
                evidence_scope=[sample_id, source], interval_convention='zero_based_half_open',
                baseline_policy=baseline, supplied_variants=supplied_variants,
                events=events, evidence_sets=evidence.sets)
