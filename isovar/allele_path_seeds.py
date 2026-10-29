"""Literal-allele witnesses in the shared original-RNA path coordinates."""

from __future__ import annotations

from collections import Counter, defaultdict
from collections.abc import Iterable, Sequence
from dataclasses import dataclass
from typing import Literal, NamedTuple

import pysam

from .allele_read import AlleleRead
from .read_collector import ReadCollector
from .read_identity import segment_identity
from .sv_rna import RnaObservation, _GeometryCache, _build_segment, _record_id
from .variant_helpers import LiteralVariantLocus


AlleleState = Literal["ALT", "REF", "OTHER"]


class OriginalAlleleCall(NamedTuple):
    """A quality-approved literal allele, even when its wider RNA is unusable."""

    state: AlleleState
    allele: str
    identity: tuple[str, str, int]
    record_id: str


@dataclass(frozen=True)
class AllelePathSeed:
    """One actual allele observation on one compatible original placement path.

    Query intervals index ``observation.sequence`` in its orientation. The
    observed interval follows the original CIGAR; the canonical interval is the
    equivalent repeat-normalized allele boundary. Both are retained, never
    confused with distinct molecules. ``record_id`` identifies the unchanged
    original SAM payload that established the allele.
    """

    observation: RnaObservation
    state: AlleleState
    allele: str
    query_interval: tuple[int, int]
    canonical_query_interval: tuple[int, int]
    record_id: str


@dataclass(frozen=True)
class AllelePathSeeds:
    """Allele/path interpretations from the supplied source-scoped templates.

    ALT/REF/OTHER segment or fragment sets can overlap across inconsistent or
    alternative placements. A reverse view and a duplicate SAM record cannot
    contribute another segment. Counts describe this acquisition scope only.
    """

    seeds: tuple[AllelePathSeed, ...]
    calls: tuple[OriginalAlleleCall, ...]
    records_considered: int
    original_records: tuple[tuple[str, str], ...]
    notes: tuple[tuple[str, int], ...]

    def segments(self, state: AlleleState) -> frozenset[tuple]:
        """All called original segments, including ones without a usable path."""
        return frozenset(call.identity for call in self.calls if call.state == state)

    def seed_segments(self, state: AlleleState) -> frozenset[tuple]:
        """Called originals that also provide an eligible allele-bearing path."""
        return frozenset(seed.observation.identity for seed in self.seeds if seed.state == state)

    def fragments(self, state: AlleleState) -> frozenset[tuple]:
        return frozenset(identity[:2] for identity in self.segments(state))


def _observation_interval(observation: RnaObservation, interval: tuple[int, int]) -> tuple[int, int] | None:
    """Map a complete sequencing-oriented interval to one observation view."""
    start, end = interval
    low, high = observation.query_interval
    if not low <= start <= end <= high:
        return None
    return (high - end, high - start) if observation.reverse_complement else (start - low, end - low)


def allele_path_seeds(
        templates: Iterable[Sequence[pysam.AlignedSegment]], locus: LiteralVariantLocus, contig: str, *,
        read_collector: ReadCollector) -> AllelePathSeeds:
    """Call alleles anew and retain their compatible observed RNA geometries.

    ``templates`` must come from one source/library and include complete
    original placement alternatives. ``contig`` is already resolved against
    that BAM header. Eligibility, repeat normalization, splice exclusion and
    technical-sequence barriers all use the existing collector. No old allele
    label, SA tag, mate position or reference base substitutes for an observed
    allele. Both orientations are retained for downstream model competition.
    """
    cache = _GeometryCache()
    reasons: Counter[str] = Counter()
    seeds: dict[tuple, AllelePathSeed] = {}
    originals: dict[str, str] = {}
    calls: set[OriginalAlleleCall] = set()
    considered = 0
    for template in templates:
        groups = defaultdict(list)
        for read in template:
            considered += 1
            if not read.query_name:
                reasons["unnamed_record"] += 1
            elif read_collector.alignment_filter_reason(read) is None:
                groups[segment_identity(read)].append(read)
        for identity, members in sorted(groups.items()):
            observations = _build_segment(identity, members, read_collector, reasons, {}, cache)
            original_members = {_record_id(read): read.to_string() for read in members}
            for read in members:
                if read.reference_name != contig:
                    continue
                extracted = read_collector.locus_read_from_pysam_aligned_segment(
                    read, **locus.collection_arguments())
                if extracted is None:
                    continue
                allele = AlleleRead.from_locus_read(extracted)
                if allele is None:
                    continue
                state: AlleleState = "ALT" if allele.allele == locus.alt else "REF" if allele.allele == locus.ref else "OTHER"
                rid = _record_id(read)
                calls.add(OriginalAlleleCall(state, allele.allele, identity, rid))
                originals[rid] = original_members[rid]
                view = read_collector.read_sequence_view(read)
                # Footprints precede trimming and index original SAM SEQ.
                original_view = view._replace(start=0, end=len(view.original_sequence))
                for footprint in allele.source_allele_qualities:
                    physical = original_view.sequenced_interval(*footprint.query_interval)
                    canonical = original_view.sequenced_interval(*footprint.canonical_query_interval)
                    for observation in observations:
                        if rid not in observation.records:
                            continue
                        interval = _observation_interval(observation, physical)
                        normalized = _observation_interval(observation, canonical)
                        if interval is None or normalized is None:
                            reasons["allele_outside_usable_observation"] += 1
                            continue
                        # An empty allele needs both observed flanks, including
                        # after an internal technical barrier split.
                        if interval[0] == interval[1] and not 0 < interval[0] < len(observation.sequence):
                            reasons["allele_boundary_missing_observed_flank"] += 1
                            continue
                        seed = AllelePathSeed(observation, state, allele.allele, interval, normalized, rid)
                        seeds[observation.key, state, interval, normalized, rid] = seed
                        originals.update((record, original_members[record]) for record in observation.records)
    return AllelePathSeeds(tuple(seeds[key] for key in sorted(seeds)), tuple(sorted(calls)), considered,
                          tuple(sorted(originals.items())), tuple(sorted(reasons.items())))
