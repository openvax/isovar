"""Annotation and focal-allele comparisons anchored to a translated RNA window."""

from __future__ import annotations

from collections.abc import Sequence
from typing import TYPE_CHECKING, TypedDict

from .genetic_code import translate_cdna
from .protein_window_comparison import ProteinComparison, compare_protein_windows, reference_window_stop
from .transcript_edit_helpers import transcript_edit_from_variant

if TYPE_CHECKING:
    from .translation import Translation


class ObservedEdit(TypedDict):
    transcript_id: str
    cdna_interval: list[int]
    alt_bases: str
    origin: str
    source_event_id: str | None


class ReferenceComparison(TypedDict):
    """Annotation and same-RNA counterfactual, with separately justified origins."""

    reference_protein_window: str
    reference_stop_in_context: bool
    reference_allele_counterfactual_protein: str
    reference_allele_counterfactual_stop: bool | None
    focal_amino_acid_change: bool | None
    annotation_comparison: ProteinComparison
    focal_allele_comparison: ProteinComparison
    comparison: str
    counterfactual_scope: str
    reference_allele_length: int
    observed_allele_length: int
    amino_acid_origins: dict[str, str]
    provenance_scope: str


def reference_comparison(translation: Translation, edits: Sequence[ObservedEdit]) -> ReferenceComparison:
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
    counter_stop_in_window = reference_window_stop(counter_amino, counter_stop, limit,
                                                   translation.ends_with_stop_codon)
    annotation_window = reference_amino[:max(0, limit)]
    annotation_stop_in_window = reference_window_stop(reference_amino, reference_stop, limit,
                                                      translation.ends_with_stop_codon)
    comparison_options = dict(allele_length_delta=delta, frameshift=translation.frameshift,
                              stop=translation.ends_with_stop_codon)
    interval = (translation.mutation_start_idx, translation.mutation_end_idx)
    annotation = compare_protein_windows(annotation_window, translation.amino_acids, interval,
                                          reference_stop=annotation_stop_in_window, **comparison_options)
    focal_comparison = compare_protein_windows(counter_window, translation.amino_acids, interval,
                                                reference_stop=counter_stop_in_window, **comparison_options)
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
    by_position: dict[int, set[str]] = {}
    for i, origin in origins:
        by_position.setdefault(i, set()).add(origin)
    return dict(reference_protein_window=reference_amino,
                reference_stop_in_context=reference_stop,
                reference_allele_counterfactual_protein=counter_window,
                reference_allele_counterfactual_stop=counter_stop_in_window,
                focal_amino_acid_change=(None if focal_comparison['status'] == 'unresolved' else
                                        focal_comparison['status'] == 'amino_acid_divergent'),
                annotation_comparison=annotation, focal_allele_comparison=focal_comparison,
                comparison='annotation_reference_cDNA_in_the_matched_frame',
                counterfactual_scope='focal allele replaced on the same observed RNA background',
                reference_allele_length=len(context.sequence_at_variant_locus),
                observed_allele_length=end - begin,
                amino_acid_origins={str(i): next(iter(values)) for i, values in by_position.items()
                                    if len(values) == 1},
                provenance_scope='codons directly touched by exactly replayed transcript edits')
