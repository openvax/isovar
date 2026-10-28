"""Indexed coverage of explicit genomic regions, without contig-name inference."""

from bisect import bisect_right
from collections import defaultdict
from collections.abc import Iterable


class GenomicIntervalIndex:
    """Merge overlapping/adjacent half-open intervals once, then query in O(log n).

    Names must already be resolved against the source header. This index describes
    acquisition or geometric coverage; it does not imply sequence or read support.
    """

    def __init__(self, intervals: Iterable[tuple[str, int, int]]) -> None:
        grouped: dict[str, list[tuple[int, int]]] = defaultdict(list)
        for contig, start, end in intervals:
            if not isinstance(contig, str) or not contig or type(start) is not int or type(end) is not int or not 0 <= start < end:
                raise ValueError("Genomic intervals require a contig and nonnegative, nonempty integer coordinates")
            grouped[contig].append((start, end))
        self._blocks: dict[str, tuple[tuple[int, ...], tuple[int, ...]]] = {}
        for contig, pieces in grouped.items():
            merged: list[tuple[int, int]] = []
            for start, end in sorted(pieces):
                if merged and start <= merged[-1][1]:
                    merged[-1] = (merged[-1][0], max(merged[-1][1], end))
                else:
                    merged.append((start, end))
            self._blocks[contig] = (tuple(start for start, _ in merged), tuple(end for _, end in merged))

    def covers(self, contig: str, start: int, end: int) -> bool:
        """Whether every base in the half-open query is in the acquired region."""
        if not 0 <= start < end:
            raise ValueError("Query interval must be nonnegative and nonempty")
        starts, ends = self._blocks.get(contig, ((), ()))
        index = bisect_right(starts, start) - 1
        return index >= 0 and ends[index] >= end

    def overlaps(self, contig: str, start: int, end: int) -> bool:
        """Whether at least one base in the half-open query is in the region."""
        if not 0 <= start < end:
            raise ValueError("Query interval must be nonnegative and nonempty")
        starts, ends = self._blocks.get(contig, ((), ()))
        index = bisect_right(ends, start)
        return index < len(starts) and starts[index] < end

    def intervals(self) -> tuple[tuple[str, int, int], ...]:
        """Return a deterministic representation for acquisition receipts."""
        return tuple((contig, start, end) for contig, (starts, ends) in sorted(self._blocks.items())
                     for start, end in zip(starts, ends))
