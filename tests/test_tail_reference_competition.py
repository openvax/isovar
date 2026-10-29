"""Candidate tails cannot erase supplied genomic evidence (#638)."""

from dataclasses import replace

import pytest
import pysam

from isovar import infer_read_ends, read_sequence_view_from_alignment
from isovar.dna import reverse_complement_dna
from isovar.read_end_inference import ReadEndInferenceCache
from isovar.read_sequence_origin import GenomicQueryEvidence


BODY = 'GCTACGTCGCTG'
SEQUENCE = BODY + 'A' * 26


@pytest.mark.parametrize('reverse', [False, True])
@pytest.mark.parametrize('competitor', ['exact', 'aligned', 'partial_exact', 'reference_break'])
def test_supplied_tail_competitors_preserve_original_sequence(reverse, competitor):
    sequence = reverse_complement_dna(SEQUENCE) if reverse else SEQUENCE
    start, end = (0, 26) if reverse else (12, 38)
    evidence = GenomicQueryEvidence(exact_intervals=((start, end),), sources=('verified-reference',))
    if competitor == 'aligned':
        evidence = replace(evidence, exact_intervals=(), aligned_intervals=((start, end),))
    elif competitor == 'partial_exact':
        evidence = replace(evidence, exact_intervals=((start + 1, end - 1),))
    elif competitor == 'reference_break':
        evidence = replace(evidence, reference_breaks=(start + 10,))
    qualities = tuple(range(len(sequence)))
    view = infer_read_ends(sequence, qualities, reverse=reverse, trim_poly_a=True, genomic_evidence=evidence)
    assert view.sequence == view.original_sequence == sequence
    assert view.qualities == qualities
    assert view.genomic_evidence == evidence
    tail, = view.annotations
    assert (tail.start, tail.end) == (start, end)
    assert tail.origin == ('reference_supported' if competitor == 'exact' else 'alignment_competitor')
    assert tail.origin_sources == ('verified-reference',)


@pytest.mark.parametrize('reverse', [False, True])
@pytest.mark.parametrize('operation', ['M', '='])
def test_aligned_poly_a_body_does_not_hide_its_unaligned_extension(reverse, operation):
    sequence = SEQUENCE
    cigar = f'20{operation}18S'
    if reverse:
        sequence, cigar = reverse_complement_dna(sequence), f'18S20{operation}'
    header = pysam.AlignmentHeader.from_references(['chr1'], [1000])
    read = pysam.AlignedSegment.fromstring(
        f'read\t{16 if reverse else 0}\tchr1\t101\t60\t{cigar}\t*\t0\t0\t{sequence}\t{"I" * len(sequence)}', header)
    original = read.to_string()
    view = read_sequence_view_from_alignment(read, trim_poly_a=True)
    assert view.sequence == (reverse_complement_dna(SEQUENCE[:20]) if reverse else SEQUENCE[:20])
    assert read.to_string() == original
    tail, = view.annotations
    assert tail.origin == 'alignment_competitor'
    assert view.genomic_evidence.sources == ('input:CIGAR',)


def test_cached_tail_matches_reassess_callable_reference_evidence():
    cache = ReadEndInferenceCache()
    evidence = GenomicQueryEvidence(exact_intervals=((12, 38),), sources=('verified-reference',))
    calls = []

    def current_reference():
        calls.append(evidence)
        return evidence

    kept = infer_read_ends(SEQUENCE, trim_poly_a=True, cache=cache, genomic_evidence=current_reference)
    evidence = GenomicQueryEvidence(limitations=('reference_unavailable',))
    trimmed = infer_read_ends(SEQUENCE, trim_poly_a=True, cache=cache, genomic_evidence=current_reference)
    assert kept.sequence == SEQUENCE and trimmed.sequence == BODY
    assert len(calls) == 2 and cache.hits == 1
    assert trimmed.annotations[0].origin == 'reference_unavailable'
    assert trimmed.genomic_evidence.limitations == ('reference_unavailable',)
    # No candidate: no reference work is needed.
    infer_read_ends(BODY, genomic_evidence=current_reference)
    assert len(calls) == 2


def test_tail_reference_bounds_are_validated_and_missing_context_is_not_negative():
    invalid = GenomicQueryEvidence(exact_intervals=((12, 39),), sources=('wrong-query',))
    with pytest.raises(ValueError, match='outside the available query'):
        infer_read_ends(SEQUENCE, genomic_evidence=invalid)
    view = infer_read_ends(SEQUENCE)
    assert view.sequence == SEQUENCE
    assert view.annotations[0].origin == 'reference_unavailable'
    assert view.genomic_evidence is None
