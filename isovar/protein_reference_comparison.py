"""Annotation and focal-allele comparisons anchored to a translated RNA window."""

from .genetic_code import translate_cdna
from .transcript_edit_helpers import transcript_edit_from_variant


def reference_comparison(translation, edits):
    """Keep the annotation comparator distinct from the same-RNA counterfactual.

    Edit provenance is attached only when replaying the transcript edits exactly
    reproduces the observed RNA. Codons touched by competing origins remain
    unclassified. These are sequence comparisons, not reference-RNA observations.
    """
    orf, context = translation.variant_orf, translation.reference_context
    offset = orf.offset_to_first_complete_codon
    begin, end = orf.variant_cdna_interval_start, orf.variant_cdna_interval_end
    reference = (orf.reference_cdna_sequence_before_variant + context.sequence_at_variant_locus
                 + orf.reference_cdna_sequence_after_variant)
    counter = orf.cdna_sequence[:begin] + context.sequence_at_variant_locus + orf.cdna_sequence[end:]
    table = context.translation_table
    reference_amino, reference_stop = translate_cdna(reference[offset:], translation_table=table)
    counter_amino, counter_stop = translate_cdna(counter[offset:], translation_table=table)
    delta = len(context.sequence_at_variant_locus) - (end - begin)
    # A capped local protein is compared at the same window extent, allowing
    # the explicit focal in-frame length difference in the reference allele.
    limit = len(translation.amino_acids) + (delta // 3 if delta % 3 == 0 else 0)
    counter_window = counter_amino[:max(0, limit)]
    counter_stop_in_window = counter_stop and len(counter_amino) <= limit
    origins = []
    for transcript in context.transcripts:
        focal = transcript_edit_from_variant(context.variant, transcript)
        left = focal.cdna_start - len(orf.reference_cdna_sequence_before_variant)
        right = left + len(reference)
        selected = sorted((e for e in edits if e['transcript_id'] == transcript.id),
                          key=lambda e: e['cdna_interval'])
        cursor, assembled, labels = left, '', []
        valid = True
        for edit in selected:
            a, b = edit['cdna_interval']
            if not left <= a <= b <= right or a < cursor:
                valid = False
                break
            prefix = reference[cursor - left:a - left]
            assembled += prefix + edit['alt_bases']
            labels += ['reference'] * len(prefix) + [edit['origin']] * len(edit['alt_bases'])
            cursor = b
        assembled += reference[cursor - left:]
        labels += ['reference'] * (right - cursor)
        if valid and assembled == orf.cdna_sequence:
            for i in range(len(translation.amino_acids)):
                touched = set(labels[offset + 3 * i:offset + 3 * (i + 1)]) - {'reference'}
                if len(touched) == 1:
                    origins.append((i, next(iter(touched))))
    by_position = {}
    for i, origin in origins:
        by_position.setdefault(i, set()).add(origin)
    return dict(reference_protein_window=reference_amino,
                reference_stop_in_context=reference_stop,
                reference_allele_counterfactual_protein=counter_window,
                reference_allele_counterfactual_stop=counter_stop_in_window,
                focal_amino_acid_change=(translation.amino_acids != counter_window or
                                        translation.ends_with_stop_codon != counter_stop_in_window),
                comparison='annotation_reference_cDNA_in_the_matched_frame',
                counterfactual_scope='focal allele replaced on the same observed RNA background',
                reference_allele_length=len(context.sequence_at_variant_locus),
                observed_allele_length=end - begin,
                amino_acid_origins={str(i): next(iter(values)) for i, values in by_position.items()
                                    if len(values) == 1},
                provenance_scope='codons directly touched by exactly replayed transcript edits')
