"""Original RNA observations and lazy indexing, independent of antigen category."""

from __future__ import annotations

from collections import Counter, OrderedDict, defaultdict
from dataclasses import dataclass
from functools import cached_property
from hashlib import sha256
from heapq import merge
from itertools import chain, combinations, islice, repeat
from types import SimpleNamespace
from typing import NamedTuple, Sequence

import numpy as np

from .cell_umi import CellUmiEvidence
from .chimeric_alignment import (
    compatible_phasing_alignments,
    source_alignment_paths_from_pysam,
)
from .dna import reverse_complement_dna
from .read_end_inference import ReadSequenceView, uncovered_intervals
from .read_identity import sam_record_id, segment_identity, source_alignments_from_pysam
from .read_lineage import ReadLineage
from .rna_path_geometry import (
    PathPositions,
    PlacementIndex,
    next_position,
    overlap_offsets,
    path_breaks,
    position_runs,
    reverse_position,
)

_RECORD_BIN = 1024  # Genomic bin width of the unbuilt-record index.


_MIN_INTRON = 21  # Shorter CIGAR N gaps are deletions (STAR's alignIntronMin default).


_MAX_SEGMENT_PATHS = 64  # Alternative supplementary groupings per segment.


@dataclass(frozen=True, eq=False)
class RnaObservation:
    """One original segment's observed path and one placement/orientation choice.

    Parameters
    ----------
    key : str
        Stable observation-choice identifier ending in + or - for reciprocal
        views. It differs between placements but is never a new read identity.
    identity : tuple[str, str, int]
        Read group, QNAME and mate bits, scoped to the owning original library.
    records : tuple[str, ...]
        sam_record_id digests of the actual contributing alignment records.
    sequence : str
        Eligible observed bases in this path orientation.
    positions : tuple
        One (contig, zero-based base, strand) or None per observed base.
    breaks : tuple
        (i, j, kind) for consecutive placed bases which are not genomic
        neighbours. Kind retains CIGAR N/D or split-record evidence.
    reverse_complement : bool
        Reverse-complement view of the processed input-read orientation;
        does not establish the biological RNA strand.
    query_interval : tuple[int, int]
        Half-open interval in the original segment's sequencing coordinates,
        unchanged between reciprocal orientation views.
    missing_qualities, secondary : bool
        Whether any contributing record lacks QUAL or is secondary.
    qualities : tuple, optional
        Per-base quality values in path orientation; None marks missing quality.
    minimum_mapping_quality : int or None, optional
        Lowest contributing MAPQ other than unavailable 255.
    mapping_quality_255 : bool, optional
        Whether any contributing record reports unavailable MAPQ 255.

    Notes
    -----
    Alternative placements and reciprocal views retain the same original
    identity. Technical barriers can split eligible query intervals; they do
    not create independent reads or invent sequence across removed bases.
    """

    key: str
    identity: tuple
    records: tuple
    sequence: str
    positions: tuple
    breaks: tuple
    reverse_complement: bool
    query_interval: tuple
    missing_qualities: bool
    secondary: bool
    qualities: tuple = ()
    # Lowest MAPQ other than 255, and whether any record reported 255.
    minimum_mapping_quality: object = None
    mapping_quality_255: bool = False

    @property
    def fragment(self):
        return self.identity[:2]

    @cached_property
    def comparison_key(self):
        return self.sequence, self.positions


def _ineligible(read, collector):
    """The ReadCollector filters, applied to each record separately."""
    return collector.alignment_filter_reason(read)


def _linked_pieces(first, second, group):
    """Reuse the reciprocal SA-path validator on two actual records of one segment."""
    sources = [source_alignments_from_pysam(r, r.query_name) for r in (first, second)]
    if sources[0] == sources[1]:
        return False
    placements = defaultdict(set)
    for read in group:
        for segment, placement in source_alignments_from_pysam(read, read.query_name):
            placements[segment].add(placement)
    observations = []
    for read, source in zip((first, second), sources):
        pairs = read.get_aligned_pairs(matches_only=True)
        if not pairs:
            return False
        q, _ = pairs[len(pairs) // 2]
        observations.append(SimpleNamespace(
            source_alignments=source,
            source_alignment_paths=source_alignment_paths_from_pysam(read, source, (q, q + 1))))
    return compatible_phasing_alignments(dict(sources[0]), *observations, placements)


def _segment_groupings(group, reasons):
    """Maximal sets of mutually linked records; alternatives are never unioned."""
    links = {(i, j) for i, j in combinations(range(len(group)), 2) if _linked_pieces(group[i], group[j], group)}
    subsets = {frozenset([i]) for i in range(len(group))}
    for i in range(len(group)):
        for subset in sorted(subsets, key=sorted):
            if i not in subset and all(tuple(sorted((i, j))) in links for j in subset):
                subsets.add(subset | {i})
        if len(subsets) > _MAX_SEGMENT_PATHS:
            reasons["segment_grouping_limit"] += 1
            break
    maximal = [s for s in subsets if not any(s < other for other in subsets)]
    return [tuple(group[i] for i in sorted(s)) for s in sorted(maximal, key=sorted)]


def _record_placements(read, view, intern):
    """Actual CIGAR placements and gap kinds, in original sequencing coordinates."""
    length = read.infer_read_length()
    hard = read.cigartuples[0][1] if read.cigartuples[0][0] == 5 else 0
    start, end = view.sequenced_interval(0, len(view.sequence))
    strand = "-" if read.is_reverse else "+"

    def original(q):
        return length - 1 - hard - q if read.is_reverse else hard + q

    placements, kinds, previous = {}, {}, None
    for q, r, n, gap, _ in _aligned_blocks(read):
        for k in range(n):
            if start <= original(q + k) < end:
                position = (read.reference_name, r + k, strand)
                placements[original(q + k)] = intern.setdefault(position, position)
        if previous is not None and gap is not None:
            kinds[tuple(sorted((original(previous), original(q))))] = _gap_kind(*gap)
        previous = q + n - 1
    return placements, kinds


def _single_record_geometry(read, view, keep_clips, intern):
    """Construct single-record geometry by CIGAR runs without per-base dictionaries."""
    blocks = list(_aligned_blocks(read))
    visible = [(max(q, view.start), min(q + n, view.end), r - q)
               for q, r, n, _, _ in blocks if max(q, view.start) < min(q + n, view.end)]
    if not visible:
        return None
    first, last = (view.start, view.end) if keep_clips else (visible[0][0], visible[-1][1])
    lo, hi = view.sequenced_interval(first - view.start, last - view.start)
    contig, strand = read.reference_name, '-' if read.is_reverse else '+'
    positions, runs = [None] * (last - first), []
    for a, b, shift in visible:
        positions[a - first:b - first] = [intern.setdefault(p, p) for p in
            zip(repeat(contig), range(a + shift, b + shift), repeat(strand))]
        start, end = (last - b, last - a) if read.is_reverse else (a - first, b - first)
        runs.append((start, end, contig, a + shift, b + shift, strand))
    if read.is_reverse:
        positions.reverse()
    runs.sort()
    merged = []
    for run in runs:
        if (merged and merged[-1][1] == run[0]
                and (merged[-1][3] == run[4] if read.is_reverse else merged[-1][4] == run[3])):
            previous = merged[-1]
            merged[-1] = (previous[0], run[1], contig, min(previous[3], run[3]), max(previous[4], run[4]), strand)
        else:
            merged.append(run)
    kinds, previous = {}, None
    for q, _, n, gap, _ in blocks:
        if previous is not None and gap is not None:
            pair = (last - 1 - q, last - 1 - previous) if read.is_reverse else (previous - first, q - first)
            kinds[pair] = _gap_kind(*gap)
        previous = q + n - 1
    positions = PathPositions(positions, merged)
    breaks = tuple((a[1] - 1, b[0], kinds.get((a[1] - 1, b[0]), 'gap')) for a, b in zip(merged, merged[1:])
                   if positions[b[0]] != next_position(positions[a[1] - 1]))
    return positions, breaks, (lo, hi)


class ObservationGeometryCache(OrderedDict):
    """Bound transient geometry reuse without sharing biological witnesses.

    Notes
    -----
    The LRU cache holds at most 8192 entries and 2,000,000 cached placement
    bases. Original sequence, quality and source identities are retained on
    each observation; equal geometry never merges independent parents.
    """

    def __init__(self):
        super().__init__()
        self.bases = 0

    def get(self, key, default=None):
        if key not in self:
            return default
        self.move_to_end(key)
        return self[key]

    def __setitem__(self, key, value):
        if key in self:
            self.bases -= self._size(self[key])
        self.bases += self._size(value)
        super().__setitem__(key, value)
        self.move_to_end(key)
        while len(self) > 8192 or self.bases > 2_000_000:
            _, previous = self.popitem(last=False)
            self.bases -= self._size(previous)

    @staticmethod
    def _size(value):
        # Mirrors retain both coordinate tuples to protect against ID reuse.
        return len(value[1]) + (len(value[0]) if len(value) == 2 else 0)


def _cache_geometry(cache, key, value):
    """Bound transient reuse independently of input size and retained witnesses."""
    if cache is not None:
        if not isinstance(cache, ObservationGeometryCache) and len(cache) >= 128:
            cache.pop(next(iter(cache)))
        cache[key] = value


class _SegmentPath(NamedTuple):
    """One contiguous original-query interval with observed geometry and QUAL."""

    sequence: str
    positions: PathPositions
    breaks: tuple[tuple[int, int, str], ...]
    query_interval: tuple[int, int]
    qualities: tuple[int | None, ...]


def _segment_path(records, collector, intern, cache=None, views=None):
    """Join one segment's linked records by original query offset.

    Returns (sequence, positions, breaks, query interval, qualities), where each
    base's quality is the lowest any record reports, or a rejection reason.
    """
    if views is None:
        views = [collector.read_sequence_view(read) for read in records]
    key = (collector.use_soft_clipped_bases, tuple(
        (read.reference_name, read.reference_start, read.is_reverse, read.cigarstring,
         view.sequence if len(records) > 1 else None,
         view.sequenced_interval(0, len(view.sequence))) for read, view in zip(records, views)))
    if len(records) == 1:
        read, view = records[0], views[0]
        start, end = view.sequenced_interval(0, len(view.sequence))
        geometry = cache.get(key) if cache is not None else None
        if geometry is None:
            built = _single_record_geometry(read, view, collector.use_soft_clipped_bases, intern)
            if built is None:
                return 'no_unique_placement'
            geometry = (None, *built)
            _cache_geometry(cache, key, geometry)
        lo, hi = geometry[3]
        sequence = (reverse_complement_dna(view.sequence) if read.is_reverse else view.sequence).upper()
        sequence = sequence[lo - start:hi - start]
        if set(sequence) - set('ACGT'):
            return 'ambiguous_bases'
        qualities = view.qualities
        if qualities is None:
            qualities = (None,) * (hi - lo)
        else:
            qualities = (qualities[::-1] if read.is_reverse else qualities)[lo - start:hi - start]
        return _SegmentPath(sequence, *geometry[1:], tuple(qualities))
    if cache is not None and key in cache:
        geometry = cache[key]
        lo, hi = geometry[3]
        # Only geometry is shared. Quality gates still see this read's original
        # qualities (including unknown values), never a representative's QUAL.
        qualities = [None] * (hi - lo)
        for read, view in zip(records, views):
            if view.qualities is None:
                continue
            start, _ = view.sequenced_interval(0, len(view.sequence))
            values = view.qualities[::-1] if read.is_reverse else view.qualities
            for q in range(max(lo, start), min(hi, start + len(values))):
                value = values[q - start]
                previous = qualities[q - lo]
                qualities[q - lo] = value if previous is None else min(previous, value)
        return _SegmentPath(*geometry, tuple(qualities))
    bases, placements, sources, kinds = {}, defaultdict(set), defaultdict(set), {}
    quality_votes = defaultdict(list)
    for index, (read, view) in enumerate(zip(records, views)):
        start, _ = view.sequenced_interval(0, len(view.sequence))
        sequence = (reverse_complement_dna(view.sequence) if read.is_reverse else view.sequence).upper()
        for q, base in enumerate(sequence, start):
            if bases.setdefault(q, base) != base:
                return "conflicting_segment_sequence"
        if view.qualities is not None:
            qualities = view.qualities[::-1] if read.is_reverse else view.qualities
            for q, quality in enumerate(qualities, start):
                quality_votes[q].append(quality)
        mapping, gaps = _record_placements(read, view, intern)
        for q, position in mapping.items():
            placements[q].add(position)
            sources[q].add(index)
        kinds.update(gaps)
    # Bases two pieces place differently (junction homology) have no single
    # observed placement; they stay in the sequence, unplaced.
    unique = {q: next(iter(p)) for q, p in placements.items() if len(p) == 1}
    if not unique:
        return "no_unique_placement"
    lo, hi = (min(bases), max(bases) + 1) if collector.use_soft_clipped_bases else (min(unique), max(unique) + 1)
    if any(q not in bases for q in range(lo, hi)):
        return "unavailable_query_bases"
    sequence = "".join(bases[q] for q in range(lo, hi))
    if set(sequence) - set("ACGT"):
        return "ambiguous_bases"
    positions = PathPositions(unique.get(q) for q in range(lo, hi))
    breaks = tuple(
        (i, j, kinds.get((i + lo, j + lo), "split" if sources[i + lo].isdisjoint(sources[j + lo]) else "gap"))
        for i, j in path_breaks(positions))
    qualities = tuple(min(quality_votes[q]) if quality_votes[q] else None for q in range(lo, hi))
    _cache_geometry(cache, key, (sequence, positions, breaks, (lo, hi)))
    return _SegmentPath(sequence, positions, breaks, (lo, hi), qualities)


def _eligible_segment_parts(built: _SegmentPath, views: Sequence[ReadSequenceView]) -> list[_SegmentPath]:
    """Slice observed geometry at technical barriers without inventing molecules."""
    sequence, positions, breaks, (lo, hi), qualities = built
    barriers = []
    for view in views:
        for start, end in view.barrier_intervals:
            start, end = max(start, view.start), min(end, view.end)
            if start < end:
                barriers.append(view.sequenced_interval(start - view.start, end - view.start))
    if not barriers:
        return [built]
    parts = []
    for start, end in uncovered_intervals(lo, hi, barriers):
        left, right = start - lo, end - lo
        local_positions = PathPositions(positions[left:right])
        if not any(position is not None for position in local_positions):
            continue
        local_breaks = tuple((a - left, b - left, kind) for a, b, kind in breaks
                             if left <= a < b < right)
        parts.append(_SegmentPath(sequence[left:right], local_positions, local_breaks,
                                  (start, end), qualities[left:right]))
    return parts


def _oriented_observations(part: _SegmentPath, identity, members, intern, cache, split: bool) -> list[RnaObservation]:
    """Represent both strands of one eligible interval with unchanged identity."""
    observations = []
    sequence, positions, breaks, interval, qualities = part
    mapqs = [r.mapping_quality for r in members]
    record_ids = tuple(sorted(sam_record_id(r) for r in members))
    key = sha256("".join(record_ids).encode()).hexdigest()[:24]
    if split:
        key += ":%d-%d" % interval
    common = dict(identity=identity, records=record_ids, query_interval=interval,
                  missing_qualities=any(r.query_qualities is None for r in members),
                  secondary=any(r.is_secondary for r in members),
                  minimum_mapping_quality=min((q for q in mapqs if q != 255), default=None),
                  mapping_quality_255=255 in mapqs)
    n = len(sequence)
    observations.append(RnaObservation(key + "+", sequence=sequence, positions=positions,
                                       breaks=breaks, reverse_complement=False, qualities=qualities, **common))
    mirror_key = ('mirror', id(positions), breaks)
    mirror = cache.get(mirror_key) if cache is not None else None
    if mirror is None:
        original_positions = positions
        positions = PathPositions((None if p is None else intern.setdefault(reverse_position(p), reverse_position(p))
                                for p in reversed(positions)),
                               [(n - b, n - a, contig, low, high, '-' if strand == '+' else '+')
                                for a, b, contig, low, high, strand in reversed(position_runs(original_positions))])
        # Hold the original tuple so its id cannot be recycled while cached.
        _cache_geometry(cache, mirror_key, (original_positions, positions))
    else:
        positions = mirror[1]
    sequence = reverse_complement_dna(sequence)
    observations.append(RnaObservation(
        key + "-", sequence=sequence, positions=positions, reverse_complement=True,
        qualities=qualities[::-1],
        breaks=tuple(sorted((n - 1 - j, n - 1 - i, kind) for i, j, kind in breaks)), **common))
    return observations


def build_segment_observations(identity, records, collector, reasons, intern, cache=None):
    """Build original-segment observation alternatives in both orientations.

    Parameters
    ----------
    identity : tuple[str, str, int]
        (read group, QNAME, mate bits), as returned by segment_identity.
    records : sequence of pysam.AlignedSegment
        Actual alignment records for that one segment; SA declarations alone
        do not supply sequence. Conflicting placements remain alternatives.
    collector : ReadCollector
        Explicit read/end policy, applied to original alignments.
    reasons : mutable counter
        Receives explicit rejection/limitation counts; not RNA-support counts.
    intern : mutable mapping
        Batch-local placement tuple interner; may be cleared between batches.
    cache : ObservationGeometryCache, optional
        Batch-local geometry reuse. Does not cache read metadata or count reads.

    Returns
    -------
    list of RnaObservation
        Eligible linked-record alternatives and reciprocal orientation views.
        Both views have the same original segment identity; technical barriers
        produce separate eligible intervals, never synthetic spanning reads.
    """
    observations = []
    for members in _segment_groupings(sorted(records, key=sam_record_id), reasons):
        if len(records) > 1 and len(members) == 1 and members[0].has_tag("SA"):
            reasons["SA_declared_piece_not_linked"] += 1
        views = [collector.read_sequence_view(read) for read in members]
        built = _segment_path(members, collector, intern, cache, views)
        if isinstance(built, str):
            reasons[built] += 1
            continue
        parts = _eligible_segment_parts(built, views)
        if len(parts) != 1 or (parts and parts[0][3] != built[3]):
            reasons['internal_technical_join_not_spanning_evidence'] += 1
        for part in parts:
            observations.extend(_oriented_observations(part, identity, members, intern, cache,
                                                         split=part.query_interval != built.query_interval))
    return observations


def _gap_kind(spliced, length):
    """A reference gap between aligned bases: N only for a plausible intron."""
    return "N" if spliced and length >= _MIN_INTRON else "D"


def _aligned_blocks(read):
    """(query offset, reference start, length, preceding gap, inserted bases) per aligned block.

    Query offsets index the stored SEQ. The preceding gap is ``(spliced,
    reference length)`` of the D/N operations since the previous block, or
    ``None``; inserted bases are the I operations since then.
    """
    q, r, gap, inserted = 0, read.reference_start, None, 0
    for operation, n in read.cigartuples:
        if operation in (0, 7, 8):
            yield q, r, n, gap, inserted
            q, r, gap, inserted = q + n, r + n, None, 0
        elif operation in (1, 4):
            q += n
            inserted += n if operation == 1 else 0
        elif operation in (2, 3):
            spliced, size = gap or (False, 0)
            gap, r = (spliced or operation == 3, size + n), r + n


def alignment_gaps(read, collector=None):
    """Yield observed CIGAR joins after the configured sequence policy.

    Parameters
    ----------
    read : pysam.AlignedSegment
        One actual alignment record.
    collector : ReadCollector
        Explicit read/end policy, applied to original alignments.

    Returns
    -------
    iterator of tuple
        (left_placement, right_placement, kind, inserted_bases), retaining
        only policy-visible flanks. This does not establish SV event linkage.
    """
    gaps, last, previous_query = [], None, None
    view = (collector.read_sequence_view(read) if collector is not None
            and (collector.trim_adapters or collector.trim_poly_a) else None)
    for q, r, n, gap, inserted in _aligned_blocks(read):
        eligible = view is None or (previous_query is not None and view.contains_usable_interval(previous_query, q + 1))
        if last is not None and gap is not None and eligible:
            gaps.append(((read.reference_name, last, "+"), (read.reference_name, r, "+"), _gap_kind(*gap),
                         read.query_sequence[q - inserted:q].upper() if inserted else ''))
        last, previous_query = r + n - 1, q + n - 1
    return gaps


class ObservationChoices(dict):
    """Original observation membership with cached fragment/placement summaries.

    Notes
    -----
    Keys identify observation alternatives. fragments deduplicates original
    RG/QNAME templates; choices preserves an original segment placement choice.
    Construct a fresh mapping before reading its cached summaries.
    """

    @cached_property
    def fragments(self):
        return {o.fragment for o in self.values()}

    @cached_property
    def choices(self):
        choices = {}
        for o in self.values():
            if o.identity not in choices or o.key < choices[o.identity]:
                choices[o.identity] = o.key
        return choices


def group_observation_choices(extenders):
    """Share overlap geometry while retaining every original member.

    Parameters
    ----------
    extenders : iterable of tuple
        (RnaObservation, relative_path_offset) overlap candidates.

    Returns
    -------
    list of tuple
        (representative, offset, ObservationChoices). A representative shares
        matching work only; fragment votes still come from all original members.
    """
    groups = {}
    for observation, offset in extenders:
        key = observation.comparison_key, offset
        if key not in groups:
            groups[key] = (observation, offset, ObservationChoices())
        groups[key][2][observation.key] = observation
    return list(groups.values())


class ObservationIndex(PlacementIndex[str]):
    """Index a frozen collection of original observation alternatives.

    Parameters
    ----------
    observations : iterable of RnaObservation
        Original observations, including reciprocal views needed for two-ended
        extension. Distinct original parents keep distinct keys.

    Notes
    -----
    grouped_overlapping returns geometric groups with their original members.
    The index is tied to this collection; do not mutate observation geometry.
    """

    def __init__(self, observations):
        self.observations = {o.key: o for o in observations}
        super().__init__({key: o.positions for key, o in self.observations.items()})
        self.windows = {}

    def grouped_overlapping(self, positions, lo, hi):
        window = positions[lo:hi]
        if window not in self.windows:
            self.windows[window] = group_observation_choices(
                (self.observations[key], offset) for key, offset in self.overlapping(window, 0, len(window)))
        return [(o, lo + offset, members) for o, offset, members in self.windows[window]]


class RnaObservationStore:
    """Lazily build original observations only when path extension needs them.

    Parameters
    ----------
    records : iterable of pysam.AlignedSegment
        Original alignment payloads; alternatives remain attached to segments.
    collector : ReadCollector
        Explicit read/end policy, applied to original alignments.
    cap : int
        Maximum overlapping original segments built for one extension query.
        limited records when a query reaches this discovery bound.
    header : dict
        Original SAM header for read-group and lineage interpretation.
    sample_id, source : str
        Explicit evidence scope, kept separate across original libraries.
    geometry_cache : ObservationGeometryCache, optional
        Transient geometry cache owned by the caller or this store.
    prebuilt : mapping, optional
        Segment identities to already-built observations from these records.

    Notes
    -----
    This store owns discovery routing and metadata, not SV classification.
    A discovery cap cannot be interpreted as a complete-input support count.
    """

    def __init__(self, records, collector, cap, header, sample_id, source, *, geometry_cache=None, prebuilt=None):
        self.collector, self.cap, self.limited = collector, cap, False
        self.sequence_overlaps = False
        self.groups, self.excluded, self.reasons = defaultdict(list), Counter(), Counter()
        self.extent, self.bins = defaultdict(dict), defaultdict(list)
        self.built, self.observations, self.joins, self.runs = {}, {}, defaultdict(list), {}
        self.intern, self.windows, self.geometry_runs = {}, {}, {}
        self.geometry_cache = ObservationGeometryCache() if geometry_cache is None else geometry_cache
        for read in records:
            reason = _ineligible(read, collector)
            if reason:
                self.excluded[reason] += 1
                continue
            identity = segment_identity(read)
            self.groups[identity].append(read)
            contig = read.reference_name
            for start, end in read.get_blocks():
                for b in range(start // _RECORD_BIN, (end - 1) // _RECORD_BIN + 1):
                    self.bins[contig, b].append((start, end, identity))
            low, high = self.extent[identity].get(contig, (read.reference_start, read.reference_end))
            self.extent[identity][contig] = (min(low, read.reference_start), max(high, read.reference_end))
        self.placement_groups = None
        self.lineage = ReadLineage(self.groups, header)
        self.cell_umi = CellUmiEvidence(self.groups, header, sample_id, source)
        if prebuilt is not None:
            for identity, observations in prebuilt.items():
                self._register(identity, observations)

    def _register(self, identity, observations):
        self.built[identity] = observations
        for o in observations:
            self.observations[o.key] = o
            for i, j, kind in o.breaks:
                self.joins[o.positions[i], o.positions[j]].append((o, i, j, kind))

    def build(self, identity):
        if identity not in self.built:
            self._register(identity, build_segment_observations(identity, self.groups[identity], self.collector, self.reasons,
                                                    self.intern, self.geometry_cache))
        return self.built[identity]

    def observation_runs(self, observation):
        if observation.key not in self.runs:
            if observation.positions not in self.geometry_runs:
                self.geometry_runs[observation.positions] = position_runs(observation.positions)
            self.runs[observation.key] = self.geometry_runs[observation.positions]
        return self.runs[observation.key]

    def _index_placement_groups(self):
        """Bin identical block/extent geometry once, retaining sorted original IDs."""
        blocks = defaultdict(set)
        for (contig, _), entries in self.bins.items():
            for start, end, identity in entries:
                blocks[identity].add((contig, start, end))
        grouped = {}
        self.placement_groups, intervals_by_contig = [], defaultdict(list)
        for identity, intervals in blocks.items():
            key = tuple(sorted(intervals)), tuple(sorted(self.extent[identity].items()))
            if key not in grouped:
                group_id = len(self.placement_groups)
                grouped[key] = group_id
                self.placement_groups.append(([], self.extent[identity]))
                for contig, start, end in intervals:
                    intervals_by_contig[contig].append((start, end, group_id))
            self.placement_groups[grouped[key]][0].append(identity)
        for members, _ in self.placement_groups:
            members.sort()
        self.placement_intervals = {contig: np.asarray(intervals, dtype=np.int64).T
                                    for contig, intervals in intervals_by_contig.items()}
        self.placement_counts = np.asarray([len(members) for members, _ in self.placement_groups], dtype=np.int64)
        self.placement_extents = {}
        self.placement_first_rank = np.empty(len(self.placement_groups), dtype=np.int64)
        for rank, group in enumerate(sorted(range(len(self.placement_groups)),
                                             key=lambda g: self.placement_groups[g][0][0])):
            self.placement_first_rank[group] = rank
        self.bins.clear()


    def overlapping(self, positions, lo, hi):
        """(key, path offset) of observations sharing a placed base in positions[lo:hi].

        Forked branches re-query shared windows; results are memoized.
        """
        return {(o.key, offset) for _, offset, members in self.grouped_overlapping(positions, lo, hi)
                for o in members.values()}

    def grouped_overlapping(self, positions, lo, hi):
        window = positions[lo:hi]
        if window not in self.windows:
            self.windows[window] = group_observation_choices(
                (self.observations[key], offset) for key, offset in self._overlapping(window))
        return [(o, lo + offset, members) for o, offset, members in self.windows[window]]

    def _overlapping(self, window):
        if self.placement_groups is None:
            self._index_placement_groups()
        runs = position_runs(window)
        coverage = np.zeros(len(self.placement_groups), dtype=np.int64)
        for _, _, contig, low, high, _ in runs:
            blocks = self.placement_intervals.get(contig)
            if blocks is not None:
                starts, ends, groups = blocks
                overlap = np.minimum(ends, high) - np.maximum(starts, low)
                hits = overlap > 0
                # Integer accumulation preserves overlapping alternate blocks;
                # this is exactly the original per-segment coverage score.
                np.add.at(coverage, groups[hits], overlap[hits])
        covered = np.flatnonzero(coverage)
        count = self.placement_counts[covered].sum()
        if count > self.cap:
            # Extenders must overlap the whole window (``min_overlap``): prefer
            # segments covering most of it (e.g. both sides of a junction),
            # then those reaching furthest beyond it.
            self.limited = True
            _, _, contig, _, _, strand = runs[-1]

            if contig not in self.placement_extents:
                infinity = np.iinfo(np.int64).max
                self.placement_extents[contig] = np.asarray(
                    [extent.get(contig, (infinity, -infinity)) for _, extent in self.placement_groups],
                    dtype=np.int64).T
            lows, highs = self.placement_extents[contig]
            reach = highs if strand == '+' else -lows
            ordered = covered[np.lexsort((-reach[covered], -coverage[covered]))]
            scores, reaches = coverage[ordered], reach[ordered]
            boundaries = np.flatnonzero((scores[1:] != scores[:-1]) | (reaches[1:] != reaches[:-1])) + 1
            candidates = []
            start = 0
            for end in chain(boundaries, [len(ordered)]):
                tied, start = ordered[start:end], end
                remaining = self.cap - len(candidates)
                if len(tied) > remaining:
                    # A group whose smallest ID trails `remaining` other group
                    # minima cannot contribute to the first `remaining` IDs.
                    chosen = np.argpartition(self.placement_first_rank[tied], remaining - 1)[:remaining]
                    tied = tied[chosen]
                candidates.extend(islice(merge(*(self.placement_groups[g][0] for g in tied)), remaining))
                if len(candidates) == self.cap:
                    break
        else:
            candidates = [identity for g in covered for identity in self.placement_groups[g][0]]
        found = set()
        for identity in candidates:
            for o in self.build(identity):
                found |= overlap_offsets(runs, 0, o.key, self.observation_runs(o))
        return found
