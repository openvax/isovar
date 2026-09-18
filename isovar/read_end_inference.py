"""Infer terminal sequence structure without modifying original alignments.

All intervals index the available original sequence, using ordinary half-open
integer coordinates. In BAM these are SAM SEQ coordinates, not forward-genomic
or hard-clip-inclusive coordinates. Trimming is one contiguous derived slice.
"""

from dataclasses import asdict, dataclass
from functools import cached_property
from hashlib import sha256
import json
import math
from pathlib import Path
from typing import NamedTuple, Optional, Tuple

from .default_parameters import (
    READ_END_WINDOW, MIN_ADAPTER_OVERLAP, MAX_ADAPTER_ERROR_RATE,
    MIN_POLY_A_LENGTH, MAX_POLY_A_ERROR_RATE, TRIM_ADAPTERS, TRIM_POLY_A,
)


def reverse_complement(sequence):
    return sequence.translate(str.maketrans("ACGTNacgtn", "TGCANtgcan"))[::-1]


@dataclass(frozen=True)
class Adapter:
    """Explicit adapter in original sequencing orientation, not BAM orientation.

    ``end`` is left or right; terminal partial suffixes/prefixes respectively
    are allowed. ``mate`` is 0 (either/unpaired), 1 or 2. Sequence must be literal
    A/C/G/T; variable barcode positions are not informative adapter sequence.
    """

    name: str
    sequence: str
    end: str = "right"
    mate: int = 0
    min_overlap: int = MIN_ADAPTER_OVERLAP
    max_error_rate: float = MAX_ADAPTER_ERROR_RATE

    def __post_init__(self):
        object.__setattr__(self, "sequence", self.sequence.upper())
        if (not self.name or self.end not in ("left", "right")
                or self.mate not in (0, 1, 2)
                or type(self.min_overlap) is not int
                or not 4 <= self.min_overlap <= len(self.sequence)
                or set(self.sequence) - set("ACGT")
                or not 0 <= self.max_error_rate < 0.5):
            raise ValueError("Invalid adapter name, sequence, end, mate or matching threshold")


@dataclass(frozen=True)
class ReadEndProfile:
    """Versioned library configuration; unknown chemistry needs no guessed kit.

    Only explicit adapter sequences are searched. ``both_orientations`` supports
    libraries with either cDNA orientation; it is not inferred from instrument.
    Poly-A/T candidates do not establish biological orientation or full length.
    """

    name: str = "unknown"
    version: str = "1"
    adapters: Tuple[Adapter, ...] = ()
    both_orientations: bool = False
    end_window: int = READ_END_WINDOW
    min_poly_a_length: int = MIN_POLY_A_LENGTH
    max_poly_a_error_rate: float = MAX_POLY_A_ERROR_RATE

    def __post_init__(self):
        object.__setattr__(self, "adapters", tuple(self.adapters))
        if (not self.name or not self.version
                or type(self.min_poly_a_length) is not int
                or type(self.end_window) is not int
                or not 4 <= self.min_poly_a_length <= self.end_window
                or not 0 <= self.max_poly_a_error_rate < 0.5
                or any(not isinstance(a, Adapter) for a in self.adapters)
                or any(len(a.sequence) > self.end_window for a in self.adapters)):
            raise ValueError("Invalid read-end profile or matching thresholds")

    @classmethod
    def from_dict(cls, values):
        values = dict(values)
        values["adapters"] = tuple(Adapter(**a) for a in values.get("adapters", ()))
        return cls(**values)

    @cached_property
    def fingerprint(self):
        """Identify the actual thresholds/sequences, not just a user-given name."""
        return sha256(json.dumps(asdict(self), sort_keys=True).encode()).hexdigest()


UNKNOWN_PROFILE = ReadEndProfile()


def read_end_profiles_from_json(filename):
    """Load one explicit profile, or ``{read_groups: {RG: profile}}``.

    Missing groups use the unknown profile; there is no implicit cross-library
    fallback. An empty-string group identifies records without an RG tag.
    """
    values = json.loads(Path(filename).expanduser().read_text())
    if "read_groups" in values:
        if set(values) != {"read_groups"}:
            raise ValueError("read_groups cannot be combined with a global profile")
        return {group: ReadEndProfile.from_dict(profile)
                for group, profile in values["read_groups"].items()}
    return ReadEndProfile.from_dict(values)


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
    original_qualities: Optional[tuple]
    start: int
    end: int
    annotations: tuple = ()
    profile: tuple = ("unknown", "1", UNKNOWN_PROFILE.fingerprint)
    reverse: bool = False
    hard_clips: tuple = (0, 0)
    upstream_poly_a_length: Optional[int] = None

    @property
    def sequence(self):
        return self.original_sequence[self.start:self.end]

    @property
    def qualities(self):
        if self.original_qualities is None:
            return None
        return self.original_qualities[self.start:self.end]

    def original_interval(self, start, end):
        """Map retained-query offsets back to available original SAM SEQ."""
        if not 0 <= start <= end <= self.end - self.start:
            raise ValueError("Interval lies outside the retained sequence")
        return self.start + start, self.start + end

    def sequenced_interval(self, start, end):
        """Map to original sequencing orientation, including unavailable H bases."""
        start, end = self.original_interval(start, end)
        left, right = self.hard_clips
        start, end = start + left, end + left
        if self.reverse:
            length = left + len(self.original_sequence) + right
            return length - end, length - start
        return start, end


def adapter_matches(sequence, adapter, window):
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
            matches.append(ReadEndAnnotation(start, end, "adapter", adapter.name,
                                              adapter.end, errors, overlap, partial))
    return tuple(matches)


def poly_a_matches(sequence, profile, start=0, end=None):
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


def infer_read_ends(sequence, qualities=None, profile=None, *, reverse=False,
                    mate=0, hard_clips=(0, 0), soft_clip_bounds=None,
                    trim_adapters=TRIM_ADAPTERS, trim_poly_a=TRIM_POLY_A,
                    upstream_poly_a_length=None):
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
    trim_adapters, trim_poly_a : bool
        Opt-in trimming. Poly-A/T inference is heuristic; this is not proof of a
        non-templated tail. Ambiguous adapter boundaries are not trimmed.
    """
    profile = UNKNOWN_PROFILE if profile is None else profile
    if qualities is not None and len(qualities) != len(sequence):
        raise ValueError("Sequence and quality lengths differ")
    if soft_clip_bounds is not None and not 0 <= soft_clip_bounds[0] <= soft_clip_bounds[1] <= len(sequence):
        raise ValueError("Invalid aligned query bounds")
    query = sequence.upper()
    candidates = []
    for adapter in profile.adapters:
        if adapter.mate and adapter.mate != mate:
            continue
        orientations = (False, True) if profile.both_orientations else (False,)
        for other_orientation in orientations:
            transformed = reverse != other_orientation
            spec = (Adapter(adapter.name, reverse_complement(adapter.sequence),
                            "left" if adapter.end == "right" else "right",
                            adapter.mate, adapter.min_overlap, adapter.max_error_rate)
                    if transformed else adapter)
            candidates.extend(adapter_matches(query, spec, profile.end_window))
    annotations = []
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
    annotations.extend(poly_a_matches(query, profile, body_start, body_end))
    start, end = 0, len(sequence)
    for annotation in annotations:
        if annotation.ambiguous:
            continue
        if annotation.kind == "adapter" and trim_adapters:
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
                            (profile.name, profile.version, profile.fingerprint), reverse, tuple(hard_clips),
                            upstream_poly_a_length)


def read_sequence_view_from_alignment(read, profile=None, *, trim_adapters=TRIM_ADAPTERS,
                                      trim_poly_a=TRIM_POLY_A):
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
        upstream_poly_a_length=read.get_tag("pt") if read.has_tag("pt") else None)
