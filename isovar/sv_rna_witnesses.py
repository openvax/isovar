"""Exact original-read witnesses for junction-linked RNA intervals."""

from bisect import bisect_right
from itertools import chain, islice

from .rna_path_geometry import position_runs

PLACEMENT_INDEX_CACHE_KEY = object()


SEQUENCE_MATCHES_CACHE_KEY = object()


_WITNESS_CACHE_ENTRIES = 4096


_WITNESS_CACHE_BASES = 250_000


def _sequence_offsets(sequence, dna):
    offset = sequence.find(dna)
    while offset >= 0:
        yield offset
        offset = sequence.find(dna, offset + 1)


class SequenceMatchCache:
    """Batch-local exact offsets shared by identical candidate DNA queries.

    Parameters
    ----------
    sequences : iterable[str]
        Frozen original observation sequences for one support batch.
    """

    def __init__(self, sequences):
        self.sequences = set(sequences)
        self.matches, self.entries, self.bases = {}, 0, 0

    def find(self, dna):
        if dna in self.matches:
            return self.matches[dna]
        room = _WITNESS_CACHE_ENTRIES - self.entries
        if room <= 0 or self.bases + len(dna) > _WITNESS_CACHE_BASES:
            return None
        matches, size = {}, 1
        for text in self.sequences:
            found = tuple(islice(_sequence_offsets(text, dna), room + 1))
            if found:
                size += len(found) + 1
                if size > room:
                    self.matches[dna] = None
                    self.entries += 1
                    self.bases += len(dna)
                    return None
                matches[text] = found
        self.matches[dna] = matches
        self.entries += size
        self.bases += len(dna)
        return matches


class PlacementCompatibility:
    """Compare placed coordinates by linear runs, allowing unplaced bases.

    Parameters
    ----------
    positions : sequence of Position or None
        One placement per observed query base; unknown placements remain None.
    """

    def __init__(self, positions):
        # Retain positions to prevent reuse of an identity-based cache key.
        self.positions = positions
        self.runs = position_runs(positions)
        self.ends = [run[1] for run in self.runs]

    def agrees(self, other, start, end, shift):
        """Every base placed on both sides has the same genomic placement."""
        i, j = bisect_right(self.ends, start), bisect_right(other.ends, start + shift)
        while i < len(self.runs) and j < len(other.runs):
            p0, p1, pc, pl, ph, ps = self.runs[i]
            q0, q1, qc, ql, qh, qs = other.runs[j]
            left, right = max(start, p0, q0 - shift), min(end, p1, q1 - shift)
            if left < right:
                if pc != qc or ps != qs:
                    return False
                own = pl + left - p0 if ps == '+' else ph - 1 - left + p0
                observed = ql + left + shift - q0 if qs == '+' else qh - 1 - left - shift + q0
                if own != observed:
                    return False
            if min(p1, q1 - shift) >= end:
                break
            if p1 <= q1 - shift:
                i += 1
            if q1 - shift <= p1:
                j += 1
        return True


def junction_interval_witnesses(sequence, positions, start, end, junctions, observations, group_cache=None):
    """Exact full-interval witnesses with compatible placements at this join.

    Parameters
    ----------
    sequence : str
        Observed RNA bases in path orientation.
    positions : sequence of Position or None
        One placement per observed query base; unknown placements remain None.
    start, end : int
        Half-open candidate interval requiring a whole original-read witness.
    junctions : sequence of dict
        Relevant junction intervals, placed flanks and direct observation keys.
    observations : mapping
        Observation key to original RnaObservation, including placement alternatives.
    group_cache : dict, optional
        Cache for one frozen observation collection; sequence and placement
        cache keys are declared by this module. Never reuse across changed observations.

    Returns
    -------
    list of dict
        Exact full-interval original witnesses and their junction links.

    Notes
    -----
    Sequence agreement includes unplaced bases. Two partial reads never
    become a full witness. Different assignments of junction homology are
    allowed only at unplaced bases, with shared placed anchors on each
    available side. A shorter ORF may need additional sequence to link it to
    the other flank; that interval must also match in this same observation.
    """
    # One cache belongs to one frozen observation store. Quality/metadata do
    # not affect sequence matching and are read from each original afterwards.
    group_cache = {} if group_cache is None else group_cache
    dna = sequence[start:end]
    search = group_cache.get(SEQUENCE_MATCHES_CACHE_KEY)
    shared_matches = search.find(dna) if search is not None else None
    if shared_matches == {}:
        return []
    witnesses = {}
    matches, consistent_intervals = {}, {}
    cached_matches = 0

    def offsets(text):
        nonlocal cached_matches
        if shared_matches is not None:
            return iter(shared_matches.get(text, ()))
        if text in matches:
            return iter(matches[text])
        iterator = _sequence_offsets(text, dna)
        room = _WITNESS_CACHE_ENTRIES - cached_matches
        if room <= 0:
            return iterator
        found = tuple(islice(iterator, room))
        if len(found) < room:
            matches[text] = found
            cached_matches += len(found) + 1
            return iter(found)
        # Repetitive sequences still yield every occurrence, without filling
        # an unbounded offsets table when no occurrence has valid placements.
        return chain(found, iterator)

    def index(pos):
        key = PLACEMENT_INDEX_CACHE_KEY, id(pos)
        if key not in group_cache:
            group_cache[key] = PlacementCompatibility(pos)
        return group_cache[key]

    for junction in junctions:
        # Flanking placed bases around the junction's unplaced interval.
        left, right = junction["query_interval"][0] - 1, junction["query_interval"][1]
        keys = tuple(junction["direct_observations"])
        if keys not in group_cache:
            groups = {}
            for key in keys:
                observation = observations[key]
                groups.setdefault(observation.comparison_key, []).append(observation)
            group_cache[keys] = list(groups.values())
        for members in group_cache[keys]:
            observation = members[0]
            if shared_matches is not None and observation.sequence not in shared_matches:
                continue
            comparison = observation.comparison_key
            for offset in offsets(observation.sequence):
                shift = offset - start
                core = comparison, shift, start, end
                if core in consistent_intervals:
                    compatible = consistent_intervals[core]
                else:
                    compatible = index(positions).agrees(index(observation.positions), start, end, shift)
                    if len(consistent_intervals) < _WITNESS_CACHE_ENTRIES:
                        consistent_intervals[core] = compatible
                if not compatible:
                    continue
                lo, hi = max(0, -shift), min(len(sequence), len(observation.sequence) - shift)

                def shared(q):
                    return positions[q] is not None and positions[q] == observation.positions[q + shift]

                anchors = []
                if junction["left"] is not None:
                    anchors.append(next((q for q in range(min(left, end - 1), lo - 1, -1) if shared(q)), None))
                if junction["right"] is not None:
                    anchors.append(next((q for q in range(max(right, start), hi) if shared(q)), None))
                if anchors and all(q is not None for q in anchors):
                    a, b = min(start, *anchors), max(end, *(q + 1 for q in anchors))
                    interval = comparison, shift, a, b
                    if interval in consistent_intervals:
                        consistent = consistent_intervals[interval]
                    else:
                        consistent = (
                            sequence[a:b] == observation.sequence[a + shift:b + shift]
                            and index(positions).agrees(index(observation.positions), a, b, shift))
                        if len(consistent_intervals) < _WITNESS_CACHE_ENTRIES:
                            consistent_intervals[interval] = consistent
                else:
                    consistent = False
                if consistent:
                    link_interval = [a + shift, b + shift]
                    for original in members:
                        a = (original.query_interval[1] - offset - len(dna) if original.reverse_complement
                             else original.query_interval[0] + offset)
                        witness = witnesses.setdefault((original.key, offset), dict(
                            observation=original.key, observation_interval=[offset, offset + len(dna)],
                            original_query_interval=[a, a + len(dna)], junction_links=[]))
                        link = dict(junction_query_interval=list(junction["query_interval"]),
                                    observation_interval=list(link_interval))
                        if link not in witness["junction_links"]:
                            witness["junction_links"].append(link)
    return [witnesses[key] for key in sorted(witnesses)]
