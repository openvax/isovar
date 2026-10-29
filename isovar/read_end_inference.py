"""Infer terminal sequence structure without modifying original alignments.

All intervals index the available original sequence, using ordinary half-open
integer coordinates. In BAM these are SAM SEQ coordinates, not forward-genomic
or hard-clip-inclusive coordinates. Trimming is one contiguous derived slice.
"""

from dataclasses import replace
from collections import OrderedDict
import math
from typing import NamedTuple, Sequence

import pysam

from .dna import reverse_complement_dna
from .default_parameters import TRIM_ADAPTERS, TRIM_POLY_A
from .read_end_profiles import (  # Compatibility imports for existing callers.
    Adapter, ReadEndProfile, UNKNOWN_PROFILE,
    read_end_profiles_from_json as read_end_profiles_from_json,
)


class ReadEndAnnotation(NamedTuple):
    """Candidate interval and reproducible match evidence, not a probability."""

    start: int
    end: int
    kind: str
    name: str
    side: str
    errors: int
    overlap: int
    partial: bool
    ambiguous: bool = False


class ReadSequenceView(NamedTuple):
    """Immutable original sequence plus a contiguous retained interval.

    Qualities are preserved as supplied, including entirely missing QUAL.
    Adapter/tail annotations always index ``original_sequence``, not ``sequence``.
    Hard clips describe unavailable bases and are not fabricated in either view.
    """

    original_sequence: str
    original_qualities: tuple[int, ...] | None
    start: int
    end: int
    annotations: tuple[ReadEndAnnotation, ...] = ()
    profile: tuple[str, str, str] = ("unknown", "1", UNKNOWN_PROFILE.fingerprint)
    reverse: bool = False
    hard_clips: tuple[int, int] = (0, 0)
    upstream_poly_a_length: int | None = None

    @property
    def sequence(self) -> str:
        return self.original_sequence[self.start:self.end]

    @property
    def qualities(self) -> tuple[int, ...] | None:
        if self.original_qualities is None:
            return None
        return self.original_qualities[self.start:self.end]

    def original_interval(self, start: int, end: int) -> tuple[int, int]:
        """Map retained-query offsets back to available original SAM SEQ."""
        if not 0 <= start <= end <= self.end - self.start:
            raise ValueError("Interval lies outside the retained sequence")
        return self.start + start, self.start + end

    def sequenced_interval(self, start: int, end: int) -> tuple[int, int]:
        """Map to original sequencing orientation, including unavailable H bases."""
        start, end = self.original_interval(start, end)
        left, right = self.hard_clips
        start, end = start + left, end + left
        if self.reverse:
            length = left + len(self.original_sequence) + right
            return length - end, length - start
        return start, end


def adapter_matches(sequence: str, adapter: Adapter, window: int) -> tuple[ReadEndAnnotation, ...]:
    """End-window edit-distance matches, including correctly anchored partials.

    Edlib provides the established bounded Levenshtein alignment. Left ends are
    reversed (not complemented) to reuse the right-end prefix matching logic.
    """
    import edlib

    reverse = adapter.end == "left"
    query = sequence[::-1] if reverse else sequence
    pattern = adapter.sequence[::-1] if reverse else adapter.sequence
    target = query[-window:]
    offset = len(query) - len(target)
    if not target:
        return ()
    matches = []
    for overlap in range(adapter.min_overlap, len(pattern) + 1):
        observed = pattern[:overlap]
        # A/T-rich and other low-information motif fragments are not adapters.
        if len(set(observed)) < 3:
            continue
        partial = overlap < len(pattern)
        budget = math.floor(overlap * adapter.max_error_rate)
        if partial:
            result = edlib.align(observed[::-1], target[::-1], mode="SHW",
                                 task="locations", k=budget)
        else:
            result = edlib.align(observed, target, mode="HW", task="locations", k=budget)
        errors = result["editDistance"]
        if errors < 0:
            continue
        locations = result["locations"]
        if not partial and errors:
            # HW reports all optimal ends, but only one start per end. Anchor
            # each end in reverse to recover every tied start as an SHW end.
            # No alignment can consume more than overlap + errors bases.
            locations = []
            for _, end in result["locations"]:
                first = max(0, end + 1 - overlap - errors)
                anchored = edlib.align(observed[::-1], target[first:end + 1][::-1],
                                       mode="SHW", task="locations", k=errors)
                locations.extend((end - reverse_end, end) for _, reverse_end in anchored["locations"])
        for start, end in locations:
            if start is None or end < start:
                continue
            if partial:
                start, end = len(target) - 1 - end, len(target) - start
            else:
                end += 1
            start, end = offset + start, offset + end
            if reverse:
                start, end = len(sequence) - end, len(sequence) - start
            matches.append(ReadEndAnnotation(start, end, adapter.role, adapter.name,
                                              adapter.end, errors, overlap, partial))
    return tuple(matches)


def poly_a_matches(sequence: str, profile: ReadEndProfile, start: int = 0,
                   end: int | None = None) -> tuple[ReadEndAnnotation, ...]:
    """Infer locally supported A/T runs; no biological-origin claim is made."""
    end = len(sequence) if end is None else end
    matches = []
    support_length = profile.min_poly_a_length
    local_budget = math.floor(support_length * profile.max_poly_a_error_rate)
    for side in ("left", "right"):
        part = sequence[start:min(end, start + profile.end_window)] if side == "left" else sequence[max(start, end - profile.end_window):end][::-1]
        for base in "AT":
            errors = score = local_errors = 0
            best = None
            for length, observed in enumerate(part, 1):
                errors += observed != base
                local_errors += observed != base
                if length > support_length:
                    local_errors -= part[length - support_length - 1] != base
                # Require an end seed and continuous local support. Do not
                # cross real sequence just because a distant tail is long.
                if length >= support_length and local_errors > local_budget:
                    break
                score += 1 if observed == base else -4
                if (length >= profile.min_poly_a_length
                        and errors <= length * profile.max_poly_a_error_rate
                        and (best is None or score > best[0])):
                    best = (score, length, errors)
            if best is not None:
                _, length, errors = best
                left, right = (start, start + length) if side == "left" else (end - length, end)
                # Window saturation is a lower bound, not a resolved boundary.
                saturated = length == profile.end_window and end - start > length
                matches.append(ReadEndAnnotation(left, right, "poly_a", "poly-" + base,
                                                  side, errors, length, saturated, saturated))
    return tuple(matches)


class _ReadEndStructure(NamedTuple):
    annotations: tuple[ReadEndAnnotation, ...]
    body_start: int
    body_end: int


def _infer_structure(sequence: str, profile: ReadEndProfile, reverse: bool, mate: int) -> _ReadEndStructure:
    """Match sequence once, independently of placement, QUAL and trimming policy."""
    candidates: list[ReadEndAnnotation] = []
    for adapter in profile.adapters:
        if adapter.mate and adapter.mate != mate:
            continue
        orientations = (False, True) if profile.both_orientations else (False,)
        for other_orientation in orientations:
            transformed = reverse != other_orientation
            spec = (replace(adapter, sequence=reverse_complement_dna(adapter.sequence),
                            end="left" if adapter.end == "right" else "right")
                    if transformed else adapter)
            candidates.extend(adapter_matches(sequence, spec, profile.end_window))
    annotations: list[ReadEndAnnotation] = []
    body_start, body_end = 0, len(sequence)
    for side in ("left", "right"):
        hits = [hit for hit in candidates if hit.side == side]
        if not hits:
            continue
        score = max(hit.overlap - 2 * hit.errors for hit in hits)
        best = set(hit for hit in hits if hit.overlap - 2 * hit.errors == score)
        cuts = {hit.end if side == "left" else hit.start for hit in best}
        ambiguous = len(cuts) != 1
        annotations.extend(hit._replace(ambiguous=ambiguous) for hit in sorted(best))
        if not ambiguous:
            cut, = cuts
            if side == "left":
                body_start = cut
            else:
                body_end = cut
    if body_start > body_end:
        annotations = [a._replace(ambiguous=True) for a in annotations]
        body_start, body_end = 0, len(sequence)
    annotations.extend(poly_a_matches(sequence, profile, body_start, body_end))
    return _ReadEndStructure(tuple(annotations), body_start, body_end)


class ReadEndInferenceCache:
    """Bound reusable sequence matches without retaining reads or file handles.

    Keys contain the exact sequence, full profile fingerprint, mate and orientation.
    Placement, qualities and upstream tags are applied to a fresh view afterwards.
    Reusing a match cannot merge identities or turn relatives into independent reads.
    Entries and retained sequence bases are bounded; oversized reads are not cached.
    """

    def __init__(self, max_entries: int = 1024, max_bases: int = 1_000_000) -> None:
        for value in (max_entries, max_bases):
            if type(value) is not int or value < 0:
                raise ValueError('Read-end cache limits must be nonnegative integers')
        self.max_entries = max_entries
        self.max_bases = max_bases
        self.retained_bases = 0
        self.hits = 0
        self.misses = 0
        self._entries: OrderedDict[tuple[str, str, bool, int], _ReadEndStructure] = OrderedDict()

    def structure(self, sequence: str, profile: ReadEndProfile, reverse: bool, mate: int) -> _ReadEndStructure:
        """Return immutable sequence evidence; never cache a placement-specific cut."""
        key = (sequence, profile.fingerprint, reverse, mate)
        if key in self._entries:
            self.hits += 1
            self._entries.move_to_end(key)
            return self._entries[key]
        self.misses += 1
        result = _infer_structure(sequence, profile, reverse, mate)
        if self.max_entries and len(sequence) <= self.max_bases:
            self._entries[key] = result
            self.retained_bases += len(sequence)
            while len(self._entries) > self.max_entries or self.retained_bases > self.max_bases:
                previous, _ = self._entries.popitem(last=False)
                self.retained_bases -= len(previous[0])
        return result


def infer_read_ends(sequence: str, qualities: Sequence[int] | None = None,
                    profile: ReadEndProfile | None = None, *, reverse: bool = False,
                    mate: int = 0, hard_clips: tuple[int, int] = (0, 0),
                    soft_clip_bounds: tuple[int, int] | None = None,
                    trim_adapters: bool = TRIM_ADAPTERS, trim_poly_a: bool = TRIM_POLY_A,
                    upstream_poly_a_length: int | None = None,
                    cache: ReadEndInferenceCache | None = None) -> ReadSequenceView:
    """Annotate sequence ends and optionally return a trimmed derived view.

    Parameters
    ----------
    sequence : str
        Unchanged available original sequence; SAM SEQ orientation for BAM input.
    profile : ReadEndProfile or None
        Explicit library profile. None allows tail inference, not guessed adapters.
    soft_clip_bounds : (int, int) or None
        For aligned input, original query alignment start/end. No bases between
        these bounds can be trimmed, including CIGAR insertions. None denotes an
        unaligned sequence, where terminal candidate trimming is unconstrained.
    reverse : bool
        Whether SAM SEQ is reversed relative to original sequencing orientation.
    cache : ReadEndInferenceCache or None
        Optional bounded cache of sequence matches. It never caches placement
        geometry, quality scores or trimming decisions.
    trim_adapters, trim_poly_a : bool
        Opt-in trimming. Poly-A/T inference is heuristic; this is not proof of a
        non-templated tail. Ambiguous adapter boundaries are not trimmed.

    Returns
    -------
    ReadSequenceView
        Unchanged available query and qualities, candidate annotations and one
        contiguous derived interval. Hard-clipped sequence remains unavailable;
        a candidate tail does not establish translation or transcript completeness.
    """
    profile = UNKNOWN_PROFILE if profile is None else profile
    if not isinstance(sequence, str) or not isinstance(profile, ReadEndProfile):
        raise ValueError('Read-end inference requires a sequence string and ReadEndProfile')
    if type(mate) is not int or mate not in (0, 1, 2) or type(reverse) is not bool:
        raise ValueError('Read-end inference requires mate 0/1/2 and boolean orientation')
    if len(hard_clips) != 2 or any(type(value) is not int or value < 0 for value in hard_clips):
        raise ValueError('Hard clips must be two nonnegative unavailable-base counts')
    if qualities is not None and len(qualities) != len(sequence):
        raise ValueError("Sequence and quality lengths differ")
    if soft_clip_bounds is not None and not 0 <= soft_clip_bounds[0] <= soft_clip_bounds[1] <= len(sequence):
        raise ValueError("Invalid aligned query bounds")
    structure = (_infer_structure(sequence.upper(), profile, reverse, mate) if cache is None else
                 cache.structure(sequence.upper(), profile, reverse, mate))
    annotations, body_start, body_end = structure
    start, end = 0, len(sequence)
    for annotation in annotations:
        if annotation.ambiguous:
            continue
        if annotation.kind in ("adapter", "primer", "linker") and trim_adapters:
            cut = annotation.end if annotation.side == "left" else annotation.start
            # An adapter match crossing into aligned evidence is a conflict,
            # not permission to trim that evidence or rewrite its CIGAR.
            if soft_clip_bounds is not None and (
                    annotation.end > soft_clip_bounds[0] if annotation.side == "left"
                    else annotation.start < soft_clip_bounds[1]):
                continue
        elif annotation.kind == "poly_a" and trim_poly_a:
            # A tail behind an untrimmed adapter cannot be removed as an
            # internal interval: that would stitch unrelated sequence together.
            if annotation.side == "left" and start < body_start:
                continue
            if annotation.side == "right" and end > body_end:
                continue
            cut = annotation.end if annotation.side == "left" else annotation.start
            if soft_clip_bounds is not None:
                cut = min(cut, soft_clip_bounds[0]) if annotation.side == "left" else max(cut, soft_clip_bounds[1])
        else:
            continue
        if annotation.side == "left":
            start = max(start, cut)
        else:
            end = min(end, cut)
    # Entirely A/T sequences can have overlapping terminal candidates.
    end = max(start, end)
    return ReadSequenceView(sequence, None if qualities is None else tuple(qualities),
                            start, end, tuple(sorted(annotations)),
                            (profile.name, profile.version, profile.fingerprint), reverse, (hard_clips[0], hard_clips[1]),
                            upstream_poly_a_length)


def read_sequence_view_from_alignment(read: pysam.AlignedSegment, profile: ReadEndProfile | None = None, *,
                                      trim_adapters: bool = TRIM_ADAPTERS,
                                      trim_poly_a: bool = TRIM_POLY_A,
                                      cache: ReadEndInferenceCache | None = None) -> ReadSequenceView:
    """Build an original-query view without altering the pysam record or tags."""
    cigar = read.cigartuples or ()
    hard_clips = (cigar[0][1] if cigar and cigar[0][0] == 5 else 0,
                  cigar[-1][1] if cigar and cigar[-1][0] == 5 else 0)
    if read.query_sequence is None:
        raise ValueError("Cannot annotate unavailable SEQ")
    return infer_read_ends(
        read.query_sequence, read.query_qualities, profile,
        reverse=read.is_reverse, mate=1 if read.is_read1 else 2 if read.is_read2 else 0,
        hard_clips=hard_clips,
        soft_clip_bounds=None if read.is_unmapped else (read.query_alignment_start, read.query_alignment_end),
        trim_adapters=trim_adapters, trim_poly_a=trim_poly_a,
        upstream_poly_a_length=read.get_tag("pt") if read.has_tag("pt") else None, cache=cache)
