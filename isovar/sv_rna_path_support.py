"""Exact original-read support for a proposed SV RNA sequence path.

Read ends, their alternative alignments, and their mates remain separate.
An exact overlap chain is a sequence hypothesis; a mate placement never fills
unsequenced bases or supplies a reading frame.
"""

from collections import defaultdict
from dataclasses import dataclass


_COMPLEMENT = str.maketrans("ACGT", "TGCA")


def reverse_complement(sequence):
    """Return the opposite DNA strand of an uppercase read sequence."""
    return sequence.translate(_COMPLEMENT)[::-1]


@dataclass(frozen=True)
class PathMatch:
    """One maximal exact read-to-path interval, in path coordinates."""

    spot: str
    mate: int
    start: int
    end: int
    clipped: bool = False


def _path_index(path, seed_size):
    index = defaultdict(list)
    for offset in range(len(path) - seed_size + 1):
        index[path[offset:offset + seed_size]].append(offset)
    return index


def exact_path_matches(path, reads, min_match=30, seed_size=20):
    """Find maximal exact sequence matches, deduplicating mapping records.

    ``reads`` yields ``(spot, mate, sequence, aligned_start, aligned_end)``.
    The aligned interval is in the supplied sequence's orientation; use
    ``None`` for unaligned original FASTQ ends. Multiple placements of one
    end can add clipping metadata but cannot add independent read support.
    """
    if not 1 <= seed_size <= min_match:
        raise ValueError("Expected 1 <= seed_size <= min_match")
    index = _path_index(path, seed_size)
    matches = {}
    for spot, mate, sequence, aligned_start, aligned_end in reads:
        sequence = sequence.upper()
        for reverse, query in ((False, sequence), (True, reverse_complement(sequence))):
            covered = defaultdict(list)
            for query_start in range(len(query) - seed_size + 1):
                for path_start in index.get(query[query_start:query_start + seed_size], ()):
                    diagonal = path_start - query_start
                    if any(lo <= query_start and query_start + seed_size <= hi
                           for lo, hi in covered[diagonal]):
                        continue
                    left = max(0, -diagonal)
                    right = min(len(query), len(path) - diagonal)
                    # A matching seed need not imply that the entire diagonal
                    # matches: extend just this seed in both directions.
                    lo, hi = query_start, query_start + seed_size
                    while lo > left and query[lo - 1] == path[lo - 1 + diagonal]:
                        lo -= 1
                    while hi < right and query[hi] == path[hi + diagonal]:
                        hi += 1
                    covered[diagonal].append((lo, hi))
                    if hi - lo < min_match:
                        continue
                    original_lo, original_hi = ((len(query) - hi, len(query) - lo)
                                                if reverse else (lo, hi))
                    clipped = (aligned_start is not None and
                               (original_lo < aligned_start or original_hi > aligned_end))
                    key = spot, mate, lo + diagonal, hi + diagonal
                    matches[key] = matches.get(key, False) or clipped
    return sorted((PathMatch(*key, clipped=clipped) for key, clipped in matches.items()),
                  key=lambda match: (match.start, match.end, match.spot, match.mate))


def interval_spots(matches, start, end):
    """Distinct original spots with one read end covering the entire interval."""
    return sorted({match.spot for match in matches
                   if match.start <= start and match.end >= end})


def overlap_chain_search(matches, start, end, min_overlap=30, max_reads=5, max_states=10000):
    """Find a distinct-fragment tiling, retaining search-budget uncertainty.

    Search is over interval/fragment histories, not globally visited reads.
    Equivalent intervals need at most ``max_reads`` distinct fragment choices:
    any chain can exclude at most ``max_reads - 1`` of them. A coverage-only
    lower bound rejects gaps cheaply before considering fragment histories.
    The returned state distinguishes no chain within the read bound from a
    state-budget cutoff. Neither establishes biological sequence absence.
    """
    if (not 0 <= start < end or min_overlap < 1 or max_reads < 1
            or type(max_states) is not int or max_states < 1):
        raise ValueError("Invalid overlap-chain bounds")
    grouped = defaultdict(dict)
    for match in matches:
        if match.end > start and match.start < end:
            grouped[match.start, match.end].setdefault(match.spot, match)
    intervals = sorted(grouped)
    result = dict(chain=[], status="no_chain_within_read_limit", states_examined=0,
                  max_states=max_states, max_reads=max_reads, min_overlap=min_overlap)
    # Greedy reach is exact when colors/fragment identities are ignored. If
    # even that relaxation fails, no distinct-fragment chain can succeed.
    reach, cursor, best = start, 0, start
    for step in range(max_reads):
        boundary = start if step == 0 else reach - min_overlap
        while cursor < len(intervals) and intervals[cursor][0] <= boundary:
            best = max(best, intervals[cursor][1])
            cursor += 1
        if best <= reach:
            return result
        reach = best
        if reach >= end:
            break
    else:
        return result
    choices = [(a, b, sorted(grouped[a, b].values(), key=lambda m: (repr(m.spot), m.mate))[:max_reads])
               for a, b in sorted(intervals, key=lambda interval: (-interval[1], interval[0]))]
    def extensions(reach, chain, used):
        boundary = reach - min_overlap if chain else start
        for a, b, rows in choices:
            if a <= boundary and b > reach:
                for match in rows:
                    if match.spot not in used:
                        yield match

    # Lazy depth-first frames avoid materializing every branch at each level.
    stack, seen = [(start, (), frozenset(), None)], set()
    while stack:
        reach, chain, used, pending = stack[-1]
        if pending is None:
            if reach >= end:
                result.update(status="found", chain=[dict(spot=m.spot, mate=m.mate,
                                                           path_interval=[m.start, m.end]) for m in chain])
                return result
            state = reach, used
            if state in seen or len(chain) >= max_reads:
                stack.pop()
                continue
            if len(seen) >= max_states:
                result.update(status="search_limit", states_examined=len(seen))
                return result
            seen.add(state)
            result["states_examined"] = len(seen)
            pending = extensions(reach, chain, used)
            stack[-1] = reach, chain, used, pending
        match = next(pending, None)
        if match is None:
            stack.pop()
        else:
            stack.append((match.end, (*chain, match), used | {match.spot}, None))
    return result


def overlap_chain(matches, start, end, min_overlap=30, max_reads=5):
    """Return a bounded distinct-fragment sequence tiling, if found.

    Use :func:`overlap_chain_search` when distinguishing a search cutoff from
    an exhaustive failure within ``max_reads`` matters. This convenience API
    preserves the historical list result; an empty list never proves absence.
    """
    return overlap_chain_search(matches, start, end, min_overlap, max_reads)["chain"]
