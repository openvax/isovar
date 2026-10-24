"""Exact original-read support for a proposed SV RNA sequence path.

Read ends, their alternative alignments, and their mates remain separate.
An exact overlap chain is a sequence hypothesis; a mate placement never fills
unsequenced bases or supplies a reading frame.
"""

from collections import defaultdict, deque
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


def overlap_chain(matches, start, end, min_overlap=30, max_reads=5):
    """Find a distinct-spot tiling from start through end, if one exists.

    Every step must overlap the previous assembled prefix by ``min_overlap``
    exact bases. This does not establish a mature isoform or translation.
    """
    if not 0 <= start < end or min_overlap < 1 or max_reads < 1:
        raise ValueError("Invalid overlap-chain bounds")
    relevant = [match for match in matches if match.end > start and match.start < end]
    first = [match for match in relevant if match.start <= start and match.end > start]
    queue = deque((match, (match,)) for match in first)
    seen = set(first)
    while queue:
        current, chain = queue.popleft()
        reach = current.end
        if reach >= end:
            return [dict(spot=match.spot, mate=match.mate,
                         path_interval=[match.start, match.end]) for match in chain]
        if len(chain) >= max_reads:
            continue
        used = {match.spot for match in chain}
        for match in relevant:
            if (match not in seen and match.spot not in used
                    and match.start <= reach - min_overlap and match.end > reach):
                seen.add(match)
                queue.append((match, (*chain, match)))
    return []
