"""Shared zero-based RNA placements, overlap indexes and exact path comparisons."""

from __future__ import annotations

from bisect import bisect_left, bisect_right
from collections import defaultdict
from collections.abc import Hashable, Mapping
from typing import Generic, Sequence, TypeVar

from .dna import reverse_complement_dna

Position = tuple[str, int, str]


def next_position(position: Position) -> Position:
    """Return the next genomic base in this oriented path.

    Parameters
    ----------
    position : Position
        Placed base in path orientation.

    Returns
    -------
    Position
        Adjacent genomic base in the same contig/strand.
    """
    contig, base, strand = position
    return contig, base + (1 if strand == "+" else -1), strand

def reverse_position(position: Position) -> Position:
    """Reverse the strand of one placed base without moving it.

    Parameters
    ----------
    position : Position
        Placed base in path orientation.

    Returns
    -------
    Position
        Same genomic base with the opposite strand.
    """
    contig, base, strand = position
    return contig, base, "-" if strand == "+" else "+"


_BIN = 64  # Genomic bin width of the observation index.


_PlacementKey = TypeVar('_PlacementKey', bound=Hashable)


class PathPositions(tuple):
    """Immutable path placements with reusable hashing and collinear runs.

    Parameters
    ----------
    values : iterable of Position or None
        Placements in query order; materialized as a tuple.
    runs : list of tuple, optional
        Precomputed position_runs(values). Supply only for this exact geometry.

    Notes
    -----
    Cached runs describe geometry only and never substitute for an original read identity.
    """

    def __new__(cls, values, runs=None):
        result = super().__new__(cls, values)
        result.hash_value = tuple.__hash__(result)
        result.runs = runs
        return result

    def __hash__(self):
        return self.hash_value


def reverse_path(sequence, positions):
    """Reverse-complement observed RNA and its placements.

    Parameters
    ----------
    sequence : str
        Observed RNA bases in the orientation of positions; no reference filling.
    positions : sequence of Position or None
        One placement per query base: (contig, zero-based base, strand), or
        None for an observed base whose placement is unresolved.

    Returns
    -------
    tuple[str, tuple]
        Reverse-complemented bases and reversed placements with strand flipped.
        Unknown placements remain None.
    """
    return reverse_complement_dna(sequence), tuple(None if p is None else reverse_position(p) for p in reversed(positions))


def position_runs(positions):
    """Encode placed query bases as maximal collinear runs.

    Parameters
    ----------
    positions : sequence of Position or None
        One placement per query base: (contig, zero-based base, strand), or
        None for an observed base whose placement is unresolved.

    Returns
    -------
    list of tuple
        (query_start, query_end, contig, genomic_low, genomic_high, strand).
        Both intervals are half-open; genomic low/high are ordered on either strand.
        Unplaced query bases split runs and are not assigned reference coordinates.
    """
    if isinstance(positions, PathPositions) and positions.runs is not None:
        return positions.runs
    runs = []
    for q, position in enumerate(positions):
        if position is None:
            continue
        if runs and runs[-1][1] == q and position == next_position(positions[q - 1]):
            start, _, contig, low, high, strand = runs[-1]
            runs[-1] = (start, q + 1, contig, min(low, position[1]), max(high, position[1] + 1), strand)
        else:
            runs.append((q, q + 1, position[0], position[1], position[1] + 1, position[2]))
    if isinstance(positions, PathPositions):
        positions.runs = runs
    return runs


def path_breaks(positions):
    """Locate discontinuities between consecutive placed query bases.

    Parameters
    ----------
    positions : sequence of Position or None
        One placement per query base: (contig, zero-based base, strand), or
        None for an observed base whose placement is unresolved.

    Returns
    -------
    list of tuple[int, int]
        Indices of the two flanking placed bases, not half-open intervals.
        Intervening unplaced bases remain between those indices.
    """
    placed = [q for q, p in enumerate(positions) if p is not None]
    return [(i, j) for i, j in zip(placed, placed[1:]) if positions[j] != next_position(positions[i])]


def merge_intervals(intervals):
    """Merge touching or overlapping intervals without filling gaps.

    Parameters
    ----------
    intervals : iterable of tuple[int, int]
        Half-open intervals in one coordinate system.

    Returns
    -------
    list of tuple[int, int]
        Sorted, disjoint half-open intervals.
    """
    merged = []
    for start, end in sorted(intervals):
        if merged and start <= merged[-1][1]:
            merged[-1] = (merged[-1][0], max(end, merged[-1][1]))
        else:
            merged.append((start, end))
    return merged


def interval_contains(intervals, position):
    """Test membership in sorted, disjoint half-open intervals.

    Parameters
    ----------
    intervals : sequence of tuple[int, int]
        Output of merge_intervals, or equivalent ordered intervals.
    position : int
        Zero-based coordinate to test.

    Returns
    -------
    bool
        Whether an interval contains the coordinate.
    """
    i = bisect_right(intervals, (position, float("inf"))) - 1
    return i >= 0 and intervals[i][0] <= position < intervals[i][1]


def overlap_offsets(query_runs, lo, key, runs):
    """Route shared placed runs to possible relative query offsets.

    Parameters
    ----------
    query_runs : sequence of tuple
        position_runs for the query slice being searched.
    lo : int
        Start of that slice in the complete query path.
    key : hashable
        Identity of the indexed path; returned unchanged.
    runs : sequence of tuple
        position_runs for the indexed path.

    Returns
    -------
    set of tuple
        (key, offset), where offset is the indexed path start in query coordinates.
        A route shares a placed base; exact sequence agreement still needs checking.
    """
    found = set()
    for start, _, contig, low, high, strand in query_runs:
        for other, _, other_contig, other_low, other_high, other_strand in runs:
            if (other_contig, other_strand) != (contig, strand) or other_low >= high or low >= other_high:
                continue
            if strand == "+":
                shared = low if low > other_low else other_low
                path_q, observed_q = start + shared - low, other + shared - other_low
            else:
                shared = (high if high < other_high else other_high) - 1
                path_q, observed_q = start + high - 1 - shared, other + other_high - 1 - shared
            found.add((key, lo + path_q - observed_q))
    return found


class PlacementIndex(Generic[_PlacementKey]):
    """Index named path placements without merging original identities.

    Parameters
    ----------
    placements : mapping
        Hashable keys to immutable placement tuples in query order.

    Notes
    -----
    overlapping(positions, lo, hi) returns (key, offset) candidates sharing a
    placed base in the half-open query interval. It is a routing index, not a
    sequence match or unique-attribution claim.
    """

    def __init__(self, placements: Mapping[_PlacementKey, tuple[tuple[str, int, str] | None, ...]]) -> None:
        grouped = {}
        for key, positions in placements.items():
            grouped.setdefault(positions, set()).add(key)
        self.geometry = [(position_runs(positions), keys) for positions, keys in grouped.items()]
        self.bins = defaultdict(set)
        for group_id, (runs, _) in enumerate(self.geometry):
            for _, _, contig, low, high, strand in runs:
                for b in range(low // _BIN, (high - 1) // _BIN + 1):
                    self.bins[contig, strand, b].add(group_id)

    def overlapping(self, positions: Sequence[tuple[str, int, str] | None], lo: int, hi: int) -> set[tuple[_PlacementKey, int]]:
        """(key, path offset of the observation start) sharing a placed base in positions[lo:hi]."""
        runs = position_runs(positions[lo:hi])
        groups = {g for _, _, contig, low, high, strand in runs
                  for b in range(low // _BIN, (high - 1) // _BIN + 1) for g in self.bins.get((contig, strand, b), ())}
        return {(key, offset) for g in groups for _, offset in overlap_offsets(runs, lo, g, self.geometry[g][0])
                for key in self.geometry[g][1]}


class PathSpans:
    """Index co-observed path geometry for junction-context spans.

    Parameters
    ----------
    sequence : str
        Observed RNA bases in the orientation of positions; no reference filling.
    positions : sequence of Position or None
        One placement per query base: (contig, zero-based base, strand), or
        None for an observed base whose placement is unresolved.

    Notes
    -----
    Calling the index with an observation returns a half-open path interval,
    or () when intervening placement/query geometry is incompatible. Placed-base
    alleles may differ: this is a geometric span, not exact full-sequence support.
    Intervening unplaced sequence and lengths must agree. Cache lifetime is one
    frozen path and observation collection.
    """

    def __init__(self, sequence, positions):
        self.sequence, self.runs, self.placed, self.cache = sequence, defaultdict(list), [0], {}
        for q0, _, contig, low, high, strand in position_runs(positions):
            self.runs[contig, strand].append((q0, low, high))
        for position in positions:
            self.placed.append(self.placed[-1] + (position is not None))
        self.starts, self.highs = {}, {}
        for key, runs in self.runs.items():
            runs.sort(key=lambda run: run[1])
            self.starts[key] = [low for _, low, _ in runs]
            maximum, highs = -1, []
            for _, _, high in runs:
                maximum = max(maximum, high)
                highs.append(maximum)
            self.highs[key] = highs

    def overlapping_runs(self, contig, strand, low, high):
        """Find potentially overlapping runs, retaining repeats and nested spans."""
        key = contig, strand
        runs = self.runs.get(key, ())
        if not runs:
            return ()
        start = bisect_right(self.highs[key], low)
        end = bisect_left(self.starts[key], high)
        return runs[start:end]

    def __call__(self, observation):
        """Return a half-open co-observed geometric span, or an empty tuple.

        The read must observe every placed path base between them, and place
        no other base between its first and last shared ones (no skipped or
        extra exon). Intervening unplaced query spans must also agree in
        length and sequence; otherwise it contributes no interval.
        """
        if observation.key not in self.cache:
            shared, shared_runs = 0, []
            shift = None
            runs = position_runs(observation.positions)
            for own, _, contig, other_low, other_high, strand in runs:
                for q0, low, high in self.overlapping_runs(contig, strand, other_low, other_high):
                    a, b = max(low, other_low), min(high, other_high)
                    if a < b:
                        shared += b - a
                        if strand == "+":
                            start, own_start = q0 + a - low, own + a - other_low
                        else:
                            start, own_start = q0 + high - b, own + other_high - b
                        offset = own_start - start
                        if shift is not None and offset != shift:
                            # Equal intervening query lengths require one
                            # read-minus-path offset across every shared block.
                            self.cache[observation.key] = ()
                            return ()
                        shift = offset
                        shared_runs.append((start, start + b - a, own_start, own_start + b - a))
            result = ()
            if shared:
                shared_runs.sort()
                first, last = shared_runs[0][0], shared_runs[-1][1] - 1
                own_first = min(a for _, _, a, _ in shared_runs)
                own_last = max(b for _, _, _, b in shared_runs) - 1
                own_placed = sum(max(0, min(end, own_last + 1) - max(start, own_first)) for start, end, *_ in runs)
                if (shared == self.placed[last + 1] - self.placed[first] == own_placed
                        and all(start - end == own_start - own_end >= 0
                                and self.sequence[end:start] == observation.sequence[own_end:own_start]
                                for (_, end, _, own_end), (start, _, own_start, _)
                                in zip(shared_runs, shared_runs[1:]))):
                    result = (first, last + 1)
            self.cache[observation.key] = result
        return self.cache[observation.key]


def exact_overlap_agrees(sequence, positions, observation, offset, lo, hi):
    """Check observed bases and jointly placed coordinates over one overlap.

    Parameters
    ----------
    sequence : str
        Observed RNA bases in the orientation of positions; no reference filling.
    positions : sequence of Position or None
        One placement per query base: (contig, zero-based base, strand), or
        None for an observed base whose placement is unresolved.
    observation : RnaObservation
        One original segment and placement choice, retaining its identity.
    offset : int
        Observation start in path coordinates; may be negative.
    lo, hi : int
        Half-open path overlap, already clipped to both sequences.

    Returns
    -------
    bool
        Exact sequence equality and no conflicting placement. A None placement
        allows the observed base to match but does not assign it a genomic coordinate.
    """
    a, b = lo - offset, hi - offset
    if observation.sequence[a:b] != sequence[lo:hi]:
        return False
    mine, theirs = positions[lo:hi], observation.positions[a:b]
    return mine == theirs or all(x is None or y is None or x == y for x, y in zip(mine, theirs))
