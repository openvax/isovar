"""Placement-specific competitors to a technical-sequence interpretation.

Intervals index the unchanged available SAM query. CIGAR/MD evidence is a claim
from the input alignment, not an independent reference measurement. Missing MD
does not establish absence of a genomic explanation.
"""

from dataclasses import dataclass
import re
from typing import Literal, Sequence

import pysam


OriginState = Literal['unassessed', 'reference_supported', 'alignment_competitor', 'reference_unavailable']
QueryInterval = tuple[int, int]


def merge_query_intervals(intervals: Sequence[QueryInterval]) -> tuple[QueryInterval, ...]:
    """Union overlapping or touching half-open query intervals."""
    merged: list[QueryInterval] = []
    for start, end in sorted(intervals):
        if merged and start <= merged[-1][1]:
            merged[-1] = (merged[-1][0], max(merged[-1][1], end))
        else:
            merged.append((start, end))
    return tuple(merged)


@dataclass(frozen=True)
class GenomicQueryEvidence:
    """Supplied genomic competitors over original query coordinates.

    ``exact_intervals`` must cover contiguous sequence verified against a named
    reference or claimed identical by CIGAR/MD. ``aligned_intervals`` alone do
    not establish sequence identity. Sources describe how evidence was obtained;
    unavailable or malformed reference context belongs in ``limitations``.
    A caller supplying external reference evidence is responsible for matching
    its query, orientation, assembly and coordinates before constructing this.
    """

    exact_intervals: tuple[QueryInterval, ...] = ()
    aligned_intervals: tuple[QueryInterval, ...] = ()
    sources: tuple[str, ...] = ()
    limitations: tuple[str, ...] = ()
    reference_breaks: tuple[int, ...] = ()

    def __post_init__(self) -> None:
        for name in ('exact_intervals', 'aligned_intervals'):
            intervals = tuple(tuple(pair) for pair in getattr(self, name))
            if any(len(pair) != 2 or any(type(value) is not int for value in pair)
                   or not 0 <= pair[0] < pair[1] for pair in intervals):
                raise ValueError('Genomic query intervals must be nonempty, nonnegative half-open pairs')
            object.__setattr__(self, name, merge_query_intervals(intervals))
        for name in ('sources', 'limitations'):
            values = getattr(self, name)
            if isinstance(values, str) or any(not isinstance(value, str) or not value.strip() for value in values):
                raise ValueError('Genomic query evidence sources and limitations must be nonempty strings')
            object.__setattr__(self, name, tuple(values))
        boundaries = tuple(self.reference_breaks)
        if any(type(value) is not int or value < 0 for value in boundaries):
            raise ValueError('Reference breaks must be nonnegative query offsets')
        object.__setattr__(self, 'reference_breaks', tuple(sorted(set(boundaries))))
        if self.exact_intervals and not self.sources:
            raise ValueError('Exact genomic query evidence requires a source')

    def validate_length(self, length: int) -> None:
        if (any(end > length for _, end in self.exact_intervals + self.aligned_intervals)
                or any(offset > length for offset in self.reference_breaks)):
            raise ValueError('Genomic evidence lies outside the available query')

    def origin(self, start: int, end: int) -> OriginState:
        """Assess one candidate without treating unexamined sequence as negative."""
        if (any(left <= start < end <= right for left, right in self.exact_intervals)
                and not any(start < offset < end for offset in self.reference_breaks)):
            return 'reference_supported'
        if any(left < end and start < right for left, right in self.aligned_intervals):
            return 'alignment_competitor'
        return 'reference_unavailable'


def _md_is_consistent_length(read: pysam.AlignedSegment) -> bool:
    """Validate MD grammar and consumption before asking htslib to expand it."""
    md = read.get_tag('MD')
    if not isinstance(md, str) or re.fullmatch(r'[0-9]+(?:(?:[A-Z]|\^[A-Z]+)[0-9]+)*', md) is None:
        return False
    observed = sum(int(token) if token.isdigit() else len(token.removeprefix('^'))
                   for token in re.findall(r'[0-9]+|\^[A-Z]+|[A-Z]', md))
    expected = sum(length for operation, length in read.cigartuples or () if operation in (0, 2, 7, 8))
    return observed == expected and not any(operation in (6, 9) for operation, _ in read.cigartuples or ())


def genomic_query_evidence_from_alignment(read: pysam.AlignedSegment) -> GenomicQueryEvidence:
    """Interpret CIGAR equality and MD claims without changing an original record.

    M and X provide placement only; I bases remain an insertion hypothesis.
    MD-derived exact runs must agree with the observed base and cannot cross a
    genomic gap or insertion. A skipped intron is not an exact genomic interval.
    CIGAR '=' claims are retained separately if MD is unavailable or malformed.
    """
    if read.is_unmapped or not read.query_sequence:
        return GenomicQueryEvidence(limitations=('unaligned_or_sequence_unavailable',))
    aligned: list[QueryInterval] = []
    exact: list[QueryInterval] = []
    sources = ['input:CIGAR']
    query = 0
    for operation, length in read.cigartuples or ():
        if operation in (0, 1, 7, 8):
            aligned.append((query, query + length))
        if operation == 7:
            exact.append((query, query + length))
        if operation in (0, 1, 4, 7, 8):
            query += length
    if query != len(read.query_sequence):
        return GenomicQueryEvidence(limitations=('CIGAR_query_length_mismatch',))
    limitations: tuple[str, ...] = ()
    if read.has_tag('MD'):
        try:
            if not _md_is_consistent_length(read):
                raise ValueError('MD and CIGAR consumption disagree')
            pairs = read.get_aligned_pairs(with_seq=True, with_cigar=True)
            runs: list[QueryInterval] = []
            previous_query = previous_reference = None
            contradictory = False
            for offset, reference, base, operation in pairs:
                if base is not None:
                    contradictory |= (operation == 7 and base.islower()) or (operation == 8 and base.isupper())
                agrees = (offset is not None and reference is not None and base is not None
                          and operation in (0, 7) and base in 'ACGT'
                          and 0 <= offset < len(read.query_sequence)
                          and read.query_sequence[offset].upper() == base)
                if agrees and offset is not None and reference is not None:
                    if runs and previous_query == offset - 1 and previous_reference == reference - 1:
                        runs[-1] = (runs[-1][0], offset + 1)
                    else:
                        runs.append((offset, offset + 1))
                    previous_query, previous_reference = offset, reference
                else:
                    previous_query = previous_reference = None
            exact = [] if contradictory else runs
            if contradictory:
                limitations = ('CIGAR_MD_equality_claims_conflict',)
            sources.append('input:MD')
        except (ValueError, AssertionError, TypeError):
            limitations = ('MD_reference_claim_unavailable_or_malformed',)
    else:
        limitations = ('MD_reference_claim_unavailable',)
    # Reference gaps have no query bases. Keep their exact boundary explicitly
    # instead of inventing a missing query base or merging across a genomic gap.
    gaps = []
    query = 0
    for operation, length in read.cigartuples or ():
        if operation in (2, 3):
            gaps.append(query)
        if operation in (0, 1, 4, 7, 8):
            query += length
    return GenomicQueryEvidence(tuple(exact), tuple(aligned), tuple(sources), limitations, tuple(gaps))
