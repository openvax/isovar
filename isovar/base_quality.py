"""Optional reported-quality assessments, separate from raw RNA support.

No platform defaults, probability calibration, or reconstruction filtering.
Indices in an AlleleQualityFootprint refer to the original SAM SEQ/QUAL,
including for reverse alignments and repeat-normalized indels.
"""

from dataclasses import dataclass
from numbers import Integral
from typing import NamedTuple

from .read_identity import source_read_ids

STATUSES = ("passed", "failed", "unassessed")


class AlleleQualityFootprint(NamedTuple):
    """Original evidence for one alignment's normalized focal allele.

    Empty alleles have two mapped flanking anchors, with None for an absent
    anchor. Their qualities describe observed bases, not deletion confidence.
    """

    source_alignment: tuple
    allele: str
    query_interval: tuple
    canonical_query_interval: tuple
    quality_indices: tuple
    quality_scores: tuple
    role: str

    @classmethod
    def from_alignment(cls, source_alignment, sequence, qualities, reference_positions,
                       canonical_interval, original_interval=None):
        """Capture before trimming or mate merging; intervals are half-open."""
        if canonical_interval is None:
            return None
        start, end = original_interval or canonical_interval
        if not (0 <= start <= end <= len(sequence)):
            return None
        if start < end:
            indices, role = tuple(range(start, end)), "allele_bases"
        else:
            indices = tuple(i if 0 <= i < len(sequence) and reference_positions[i] is not None else None
                            for i in (start - 1, end))
            role = "flanking_anchors"
        scores = tuple(None if i is None or qualities is None else qualities[i] for i in indices)
        return cls(source_alignment, sequence[slice(*canonical_interval)], (start, end),
                   canonical_interval, indices, scores, role)


class BaseQualityAssessment(NamedTuple):
    """One source alignment's quality assessment for the observation's allele."""

    source_alignment: tuple
    status: str
    reason: str
    minimum_reported_quality: object
    missing_bases: int


@dataclass(frozen=True)
class BaseQualityPolicy:
    """Assess focal-allele reported qualities against an explicit threshold.

    Parameters
    ----------
    min_base_quality : int
        Minimum reported Phred score, inclusive, from 0 to 254. Missing QUAL
        (including BAM's 255 sentinel) is unknown. This is neither a universal
        platform cutoff nor a calibrated probability of sequence correctness.
    """

    min_base_quality: int

    def __post_init__(self):
        if (isinstance(self.min_base_quality, bool) or not isinstance(self.min_base_quality, Integral)
                or not 0 <= self.min_base_quality < 255):
            raise ValueError("min_base_quality must be an integer from 0 to 254")

    def description(self):
        """JSON-ready definition shared by all optional quality outputs."""
        return dict(min_base_quality=int(self.min_base_quality), scope="focal_allele",
                    empty_allele="mapped_flanking_anchor_quality_only", missing="unassessed",
                    conflicting_source_allele="unassessed", reconstruction_filtered=False,
                    calibrated_probability=False)

    def assess(self, read):
        """Assess original source footprints, never the merged consensus QUAL.

        Caller-created/legacy observations without original provenance remain
        unassessed. A mate with a different original allele cannot qualify as
        quality support for the merged allele.
        """
        footprints = getattr(read, "source_allele_qualities", ())
        assessments = []
        for footprint in footprints:
            known = [int(q) for q in footprint.quality_scores
                     if isinstance(q, Integral) and not isinstance(q, bool) and 0 <= q < 255]
            missing = len(footprint.quality_scores) - len(known)
            minimum = min(known) if known else None
            if footprint.allele != read.allele:
                status, reason = "unassessed", "source_allele_disagrees"
            elif minimum is not None and minimum < self.min_base_quality:
                status, reason = "failed", "below_threshold"
            elif missing or not known:
                status, reason = "unassessed", "missing_quality_or_anchor"
            else:
                status, reason = "passed", "reported_quality_meets_threshold"
            assessments.append(BaseQualityAssessment(footprint.source_alignment, status, reason, minimum, missing))
        represented = {a.source_alignment for a in assessments}
        for alignment in getattr(read, "source_alignments", ()):
            if alignment not in represented:
                assessments.append(BaseQualityAssessment(alignment, "unassessed", "missing_provenance", None, 0))
        return tuple(assessments) or (BaseQualityAssessment(None, "unassessed", "missing_provenance", None, 0),)

    def partition(self, reads):
        """Deduplicate segments conservatively across their observed placements.

        Any known failure fails; otherwise all placements must pass. Missing
        provenance stays unassessed. Returns identity sets and legacy read lists
        by status. Fragment/cell sets may overlap between these categories.
        """
        states, legacy = {}, {status: [] for status in STATUSES}
        for read in reads:
            keys = source_read_ids(read)
            if not keys:
                legacy["unassessed"].append(read)
                continue
            for assessment in self.assess(read):
                if assessment.source_alignment is not None:
                    key = assessment.source_alignment[0]
                    if key in keys:
                        states.setdefault(key, set()).add(assessment.status)
        partitions = {status: set() for status in STATUSES}
        for key, statuses in states.items():
            status = "failed" if "failed" in statuses else "passed" if statuses == {"passed"} else "unassessed"
            partitions[status].add(key)
        return partitions, legacy
