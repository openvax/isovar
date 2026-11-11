"""Original RNA coverage of an explicitly positioned coding interval."""

from bisect import bisect_left
from collections import defaultdict
from collections.abc import Callable, Iterable, Sequence
from itertools import accumulate
from typing import TypedDict

from .allele_read import AlleleRead
from .dna import reverse_complement_dna
from .read_identity import fragment_ids
from .rna_evidence import RNASupport


class CodingIntervalSupport(TypedDict):
    scope: str
    coding_interval: list[int]
    fragment_depth: list[int]
    minimum_fragment_depth: int
    fully_tiled: bool
    allele_observation_span_support: RNASupport
    single_read_span_support: RNASupport
    paired_template_span_support: RNASupport


def _union(spans):
    """Merge half-open spans into sorted disjoint ones."""
    merged = []
    for start, end in sorted(spans):
        if merged and start <= merged[-1][1]:
            merged[-1][1] = max(merged[-1][1], end)
        else:
            merged.append([start, end])
    return merged


class AlleleReadTiling:
    """Where allele observations lie on one candidate cDNA, and where they disagree with it.

    Placing a read depends only on the cDNA, the focal allele offset and the
    strand, so one tiling answers `interval_support` for every coding interval
    of that cDNA without comparing the reads again.

    Parameters
    ----------
    reads : iterable of AlleleRead
        Original allele observations, possibly merged overlapping mates.
    sequence : str
        Candidate cDNA in transcript orientation.
    variant_start : int
        Zero-based start of the focal alternate allele in ``sequence``.
    strand : str
        Transcript strand relative to the observations' genomic orientation.

    Attributes
    ----------
    sequence : str
        The candidate cDNA.
    variant_start : int
        The focal allele offset.
    """

    def __init__(self, reads: Iterable[AlleleRead], sequence: str, variant_start: int, strand: str):
        if strand not in ("+", "-"):
            raise ValueError("Invalid strand")
        if not 0 <= variant_start <= len(sequence):
            raise ValueError("Invalid focal allele offset")
        self.sequence, self.variant_start = sequence, variant_start
        # (start, end, mismatched offsets, read) of each read overlapping the cDNA, in read order.
        self._placements = []
        fragments = defaultdict(list)
        for read in reads:
            observed, allele_start = read.sequence, len(read.prefix)
            if strand == "-":
                observed, allele_start = reverse_complement_dna(observed), len(read.suffix)
            shift = variant_start - allele_start
            start, end = max(0, shift), min(len(sequence), shift + len(observed))
            if end <= start:
                continue
            mismatches = () if sequence[start:end] == observed[start - shift:end - shift] else tuple(
                offset for offset in range(start, end) if sequence[offset] != observed[offset - shift])
            placement = start, end, mismatches, read
            self._placements.append(placement)
            for fragment in fragment_ids([read]):
                fragments[fragment].append(placement)
        # A fragment whose reads all agree with the cDNA covers the same bases
        # in every interval. One with a mismatched read is admitted per
        # interval, because a mismatch outside the interval does not exclude it.
        changes = [0] * (len(sequence) + 1)
        self._conditional = []
        for placements in fragments.values():
            if any(mismatches for _, _, mismatches, _ in placements):
                self._conditional.append(placements)
                continue
            for start, end in _union((start, end) for start, end, _, _ in placements):
                changes[start] += 1
                changes[end] -= 1
        self._depth = list(accumulate(changes[:-1]))

    def interval_support(self, interval: Sequence[int],
                         support: Callable[[list[AlleleRead]], RNASupport]) -> CodingIntervalSupport:
        """Measure the reads' exact coverage of one coding interval.

        Parameters
        ----------
        interval : pair of int
            Half-open coding interval, including an observed stop when requested.
        support : callable
            Converts a read collection to its scoped RNA support record.

        Returns
        -------
        dict
            Per-base fragment depth, exact full-interval allele-observation support
            and its single-read and merged paired-template subsets. The latter
            counts observations made from two overlapping mates, not unsequenced
            inserts. Single-read spans exclude merged observations; they do not
            recount the original mates hidden within a merged consensus. Different placements
            and overlapping observations of a fragment count once at each base.
        """
        start, end = interval
        if not 0 <= start <= end <= len(self.sequence):
            raise ValueError("Invalid coding interval")

        def clipped(placement):
            """The placement's part of the interval, or None if absent or mismatched there."""
            left, right = max(start, placement[0]), min(end, placement[1])
            mismatches = placement[2]
            if left < right and bisect_left(mismatches, left) == bisect_left(mismatches, right):
                return left, right
            return None

        depths = self._depth[start:end]
        for placements in self._conditional:
            for left, right in _union(filter(None, map(clipped, placements))):
                for offset in range(left, right):
                    depths[offset - start] += 1
        spanning = [placement[3] for placement in self._placements if clipped(placement) == (start, end)]
        return dict(
            scope="exact_observed_coding_interval",
            coding_interval=[start, end], fragment_depth=depths,
            minimum_fragment_depth=min(depths) if depths else 0,
            fully_tiled=bool(depths) and all(depths),
            allele_observation_span_support=support(spanning),
            single_read_span_support=support([read for read in spanning if read.source_read_count == 1]),
            paired_template_span_support=support([read for read in spanning if read.source_read_count == 2]))


def coding_interval_support(
        reads: Iterable[AlleleRead], sequence: str, variant_start: int,
        interval: Sequence[int], strand: str,
        support: Callable[[list[AlleleRead]], RNASupport]) -> CodingIntervalSupport:
    """Measure exact allele-anchored coverage without counting partial contributors as spans.

    Parameters
    ----------
    reads : iterable of AlleleRead
        Original allele observations, possibly merged overlapping mates.
    sequence : str
        Candidate cDNA in transcript orientation.
    variant_start : int
        Zero-based start of the focal alternate allele in ``sequence``.
    interval : pair of int
        Half-open coding interval, including an observed stop when requested.
    strand : str
        Transcript strand relative to the observations' genomic orientation.
    support : callable
        Converts a read collection to its scoped RNA support record.

    Returns
    -------
    dict
        See `AlleleReadTiling.interval_support`, which several intervals of one
        cDNA should share.
    """
    return AlleleReadTiling(reads, sequence, variant_start, strand).interval_support(interval, support)
