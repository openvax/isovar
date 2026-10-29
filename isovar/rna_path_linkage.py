"""Feature linkage and complete-path coverage from original placement choices.

A match is one original segment on one alternative geometry. Never union two
placements of that segment. Paired links join different read ends of one RG /
QNAME template; they do not supply sequence between those ends.
"""

from __future__ import annotations

from collections.abc import Sequence
from dataclasses import dataclass
from itertools import combinations


@dataclass(frozen=True)
class PathFeature:
    """A query interval requiring observation, optionally by particular originals.

    Empty intervals are deletion/junction boundaries and need bases on both
    sides. ``observation_keys`` restricts a focal allele to original placements
    that passed the allele caller, including its quality policy. None permits
    any matching original; an empty tuple permits none.
    """

    id: str
    start: int
    end: int
    observation_keys: tuple[str, ...] | None = None


@dataclass(frozen=True)
class OriginalPathMatch:
    """One coordinate-consistent overlap of one original placement with a path."""

    identity: tuple[str, str, int]
    observation_key: str
    start: int
    end: int

    def covers(self, feature: PathFeature) -> bool:
        """Whether this single placement observes the complete feature."""
        allowed = feature.observation_keys
        if allowed is not None and self.observation_key not in allowed:
            return False
        if feature.start == feature.end:
            return self.start < feature.start < self.end
        return self.start <= feature.start < feature.end <= self.end


def paired_ends(first: OriginalPathMatch, second: OriginalPathMatch) -> bool:
    """Two different physical read ends of exactly one source-scoped template."""
    return first.identity[:2] == second.identity[:2] and {first.identity[2], second.identity[2]} == {64, 128}


def template_path_linkage(matches: Sequence[OriginalPathMatch], length: int,
                          features: Sequence[PathFeature]) -> dict:
    """Summarize one template without combining its alternative placements.

    Each pairwise feature link requires one compatible choice for each of at
    most two read ends. Links for different feature pairs need not share a
    globally consistent placement and cannot be combined into whole-path phase.
    Complete paired coverage requires the observed intervals to touch/overlap.
    """
    whole = PathFeature("whole_path", 0, length)
    singles = {m.identity for m in matches if m.covers(whole)}
    pairs = [(a, b) for a, b in combinations(matches, 2) if paired_ends(a, b)]
    paired_coverage = any(min(a.start, b.start) == 0 and max(a.end, b.end) == length
                          and max(a.start, b.start) <= min(a.end, b.end) for a, b in pairs)
    links = []
    for first, second in combinations(features, 2):
        direct = {m.identity for m in matches if m.covers(first) and m.covers(second)}
        paired = any((a.covers(first) and b.covers(second)) or (b.covers(first) and a.covers(second))
                     for a, b in pairs)
        if direct or paired:
            links.append(dict(features=[first.id, second.id], single_read_segments=sorted(direct),
                              paired_template=paired))
    return dict(whole_path_single_read_segments=sorted(singles),
                whole_path_paired_template=paired_coverage, links=links)
