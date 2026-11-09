"""Read and fragment support for an explicit local RNA feature.

The counting predicate and geometric sampling opportunities share one interval
implementation. Antigen categories do not enter either calculation. These are
observations of RNA sequence/path, not translation or molecule abundance.
"""

from __future__ import annotations

from collections import defaultdict
from dataclasses import asdict, dataclass
from itertools import combinations
import math
from typing import Any, Iterable

from .rna_evidence import content_identifier
from .rna_path_linkage import OriginalPathMatch, paired_ends


SCHEMA = "isovar.rna_feature_support.v1"
GEOMETRY_MODEL = "uniform_feature_start_opportunities.v1"
FEATURE_KINDS = frozenset(("allele", "junction", "coding_interval", "compound"))


def _integer(value, name, minimum=0):
    if type(value) is not int or value < minimum:
        raise ValueError(f"{name} must be an integer >= {minimum}")


def _intervals(values, name):
    result = tuple(tuple(v) for v in values)
    if not result:
        raise ValueError(f"{name} must not be empty")
    for pair in result:
        if len(pair) != 2:
            raise ValueError(f"{name} requires half-open interval pairs")
        start, end = pair
        _integer(start, name)
        _integer(end, name)
        if end <= start:
            raise ValueError(f"{name} intervals must be nonempty; expand junction anchors")
    return result


@dataclass(frozen=True)
class RNAFeature:
    """An explicit set of required intervals on one oriented RNA path.

    Parameters
    ----------
    identifier, path_id : str
        Stable feature and carrier-path identities.
    intervals : tuple of pairs of int
        Required half-open nucleotide intervals, including required anchors.
        A junction is represented by its left and right anchor bases.
    kind : str
        ``allele``, ``junction``, ``coding_interval`` or ``compound``.
    mode : str
        ``all`` requires every interval in one read/template; ``any`` is
        explicitly partial/event support. Separate templates never create phase.
    path_length : int or None
        A known complete RNA path length. None uses an explicitly unbounded
        local-path geometry, not an inferred cropped transcript end.
    observation_keys : tuple of str or None
        Optional original placement keys admitted by an upstream allele caller.
        An empty tuple permits no original; None imposes no additional filter.
    """

    identifier: str
    path_id: str
    intervals: tuple[tuple[int, int], ...]
    kind: str
    mode: str = "all"
    path_length: int | None = None
    observation_keys: tuple[str, ...] | None = None

    def __post_init__(self):
        if not self.identifier or not self.path_id or self.kind not in FEATURE_KINDS:
            raise ValueError("RNA feature requires identities and a supported kind")
        object.__setattr__(self, "intervals", _intervals(self.intervals, "feature"))
        if self.mode not in ("all", "any"):
            raise ValueError("RNA feature mode must be all or any")
        if self.path_length is not None:
            _integer(self.path_length, "path_length", 1)
            if max(b for _, b in self.intervals) > self.path_length:
                raise ValueError("RNA feature lies outside the complete path")
        if self.observation_keys is not None:
            object.__setattr__(self, "observation_keys", tuple(self.observation_keys))
            if any(not isinstance(k, str) or not k for k in self.observation_keys):
                raise ValueError("Original observation keys must be nonempty strings")


@dataclass(frozen=True)
class RNASamplingGeometry:
    """One possible sequenced-template geometry and its sampling weight.

    ``segments`` contains one SE interval or the two sequenced PE intervals
    relative to the template's start on the RNA path. The gap between PE
    segments is not observed. ``template_length`` is in RNA coordinates,
    never the genomic TLEN across an intron. Weights need not sum to one.
    """

    segments: tuple[tuple[int, int], ...]
    template_length: int
    weight: float = 1.0

    def __post_init__(self):
        object.__setattr__(self, "segments", _intervals(self.segments, "segments"))
        _integer(self.template_length, "template_length", 1)
        if len(self.segments) not in (1, 2) or max(b for _, b in self.segments) > self.template_length:
            raise ValueError("Geometry requires one or two segments within the template")
        if isinstance(self.weight, bool) or not math.isfinite(self.weight) or self.weight <= 0:
            raise ValueError("Geometry weight must be positive and finite")


@dataclass(frozen=True)
class RNASamplingProfile:
    """Library geometry and independently supplied whole-library denominators.

    Counts and geometry can be absent independently. ``provenance`` identifies
    the original library statistics and counting policy. An extract's record
    count is not a whole-library denominator. The uniform model does not correct
    sequence, mapping, capture or compositional bias.
    """

    provenance: str
    geometries: tuple[RNASamplingGeometry, ...] = ()
    library_reads: int | None = None
    library_fragments: int | None = None
    denominator_scope: str | None = None

    def __post_init__(self):
        if not isinstance(self.provenance, str) or not self.provenance:
            raise ValueError("RNA sampling profile requires provenance")
        object.__setattr__(self, "geometries", tuple(self.geometries))
        if any(not isinstance(g, RNASamplingGeometry) for g in self.geometries):
            raise TypeError("Sampling geometries must be RNASamplingGeometry values")
        for name in ("library_reads", "library_fragments"):
            if getattr(self, name) is not None:
                _integer(getattr(self, name), name, 1)
                if self.denominator_scope != "whole_library":
                    raise ValueError("RNA normalization requires a whole_library denominator")
        if self.denominator_scope not in (None, "whole_library"):
            raise ValueError("Only declared whole-library denominators are accepted")
        if self.library_reads is not None and self.library_fragments is not None:
            if not self.library_fragments <= self.library_reads <= 2 * self.library_fragments:
                raise ValueError("Read and fragment denominators disagree")
            if self.geometries and all(len(g.segments) == 1 for g in self.geometries):
                if self.library_reads != self.library_fragments:
                    raise ValueError("Single-end read and fragment denominators must agree")

    @classmethod
    def from_dict(cls, value):
        """Restore a JSON profile with strict field validation."""
        value = dict(value)
        value["geometries"] = tuple(RNASamplingGeometry(**g) for g in value.get("geometries", ()))
        return cls(**value)


def _merge(intervals):
    merged = []
    for start, end in sorted(intervals):
        if merged and start <= merged[-1][1]:
            merged[-1] = (merged[-1][0], max(end, merged[-1][1]))
        else:
            merged.append((start, end))
    return merged


def _accepts(observed, feature):
    observed = _merge(observed)
    covered = [any(start <= low and end >= high for start, end in observed)
               for low, high in feature.intervals]
    return all(covered) if feature.mode == "all" else any(covered)


def effective_rna_feature_lengths(feature: RNAFeature, profile: RNASamplingProfile) -> dict:
    """Return separate expected read/fragment opportunity lengths in effective bp.

    Enumerate changes in interval coverage rather than bases in the complete
    transcript. Read opportunities average over sequenced ends, matching the
    read denominator; fragment opportunities count the observed pair once.
    ``None`` means sampling geometry was not provided; zero means the declared
    feature is unobservable under this model. The method is category-independent.
    """
    if not profile.geometries:
        return dict(reads=None, fragments=None)
    read_total = fragment_total = read_weight = fragment_weight = 0.0
    for geometry in profile.geometries:
        segments = geometry.segments
        # Acceptance can change only when an observed endpoint crosses a
        # required endpoint. Include the integer point after each crossing.
        transitions = set()
        for low, high in feature.intervals:
            for start, end in segments:
                transitions.update((low - start, low - start + 1, high - end, high - end + 1))
        minimum = min(a for a, _ in feature.intervals) - max(b for _, b in segments)
        maximum = max(b for _, b in feature.intervals) - min(a for a, _ in segments) + 1
        if feature.path_length is not None:
            minimum = max(minimum, 0)
            maximum = min(maximum, feature.path_length - geometry.template_length + 1)
        if maximum <= minimum:
            starts = []
        else:
            points = sorted({minimum, maximum, *(x for x in transitions if minimum < x < maximum)})
            starts = zip(points, points[1:])
        read_opportunities = fragment_opportunities = 0
        for start, end in starts:
            observed = [(start + a, start + b) for a, b in segments]
            read_opportunities += (end - start) * sum(_accepts([segment], feature) for segment in observed)
            fragment_opportunities += (end - start) * _accepts(observed, feature)
        read_total += geometry.weight * read_opportunities
        read_weight += geometry.weight * len(segments)
        fragment_total += geometry.weight * fragment_opportunities
        fragment_weight += geometry.weight
    return dict(reads=read_total / read_weight, fragments=fragment_total / fragment_weight)


def _count_matches(feature, matches):
    by_fragment = defaultdict(list)
    for match in set(matches):
        if (not isinstance(match, OriginalPathMatch) or len(match.identity) != 3
                or match.identity[2] not in (0, 64, 128) or match.start < 0 or match.end <= match.start):
            raise ValueError("Feature support requires original path placements and segment identities")
        if feature.path_length is not None and match.end > feature.path_length:
            raise ValueError("Observed match extends beyond the complete RNA path")
        if feature.observation_keys is None or match.observation_key in feature.observation_keys:
            by_fragment[match.identity[:2]].append(match)
    reads, fragments, joint = set(), set(), set()
    for identity, rows in by_fragment.items():
        single = {row.identity for row in rows if _accepts([(row.start, row.end)], feature)}
        reads.update(single)
        if single:
            fragments.add(identity)
        elif any(paired_ends(a, b) and _accepts([(a.start, a.end), (b.start, b.end)], feature)
                 for a, b in combinations(rows, 2)):
            fragments.add(identity)
            joint.add(identity)
    return reads, fragments, joint


def measure_rna_feature(
        feature: RNAFeature, matches: Iterable[OriginalPathMatch], *,
        evidence_scope: tuple[str, str], profile: RNASamplingProfile | None = None,
        provenance: dict[str, Any] | None = None) -> dict:
    """Count original RNA observations and normalize the same local feature.

    ``matches`` are coordinate- and sequence-consistent original placements,
    as produced by Isovar's RNA path reconstruction. At most one alternative
    placement per read end can contribute to a template's complete support.
    An ``all`` compound feature may use two ends of one template, never two
    unrelated templates. Inputs must have undergone the declared eligibility
    policy; ``provenance`` records that policy and attribution limitations.

    Returned JSON uses hashed read/template identities, with raw counts and
    unavailable reasons. Sampling profile provenance and the feature are part
    of the measurement identity. A supplied regional set remains regional even
    when a verified original-library denominator is supplied separately.
    """
    if len(evidence_scope) != 2 or any(not isinstance(v, str) or not v for v in evidence_scope):
        raise ValueError("RNA evidence scope requires sample and source identities")
    reads, fragments, joint = _count_matches(feature, matches)
    scope = list(evidence_scope)
    profile = profile or RNASamplingProfile("not_supplied")
    lengths = effective_rna_feature_lengths(feature, profile)
    ids = dict(read_ids=sorted(content_identifier("segment", [scope, list(k)]) for k in reads),
               fragment_ids=sorted(content_identifier("fragment", [scope, list(k)]) for k in fragments),
               joint_only_fragment_ids=sorted(content_identifier("fragment", [scope, list(k)]) for k in joint))
    values, unavailable = {}, {}
    for unit, count, denominator in (("reads", len(reads), profile.library_reads),
                                     ("fragments", len(fragments), profile.library_fragments)):
        length = lengths[unit]
        if length == 0 and count:
            raise ValueError(f"Observed {unit} support contradicts zero effective opportunity")
        if denominator is not None and count > denominator:
            raise ValueError(f"Supporting {unit} exceed the declared library denominator")
        values[unit] = count
        values["library_" + unit] = denominator
        values["effective_length_" + unit] = length
        for suffix, scale, divisors, reason in (
                ("per_million", 1e6, (denominator,), "whole_library_denominator_unavailable"),
                ("per_effective_kb", 1e3, (length,), "sampling_geometry_unavailable" if length is None else "feature_unobservable"),
                ("per_effective_kb_per_million", 1e9, (denominator, length),
                 "whole_library_denominator_unavailable" if denominator is None else
                 "sampling_geometry_unavailable" if length is None else "feature_unobservable")):
            name = unit + "_" + suffix
            if any(v is None or v == 0 for v in divisors):
                values[name] = None
                unavailable[name] = reason
            else:
                values[name] = scale * count / math.prod(divisors)
    record = dict(schema=SCHEMA, feature=asdict(feature), evidence_scope=scope,
                  sampling_profile=asdict(profile), normalization_model=GEOMETRY_MODEL,
                  geometry_scope="complete_path" if feature.path_length is not None else "unbounded_local_path",
                  values=values, unavailable=unavailable, **ids,
                  provenance=dict(provenance or {}))
    record["measurement_id"] = content_identifier("rna_feature", record)
    return record


def validate_rna_feature_support(record: dict) -> dict:
    """Validate a saved measurement without recounting reads or trusting scalars.

    Recompute its identity, geometry, arithmetic and identity cardinalities.
    Return the original JSON-compatible record for embedding in consumer data.
    """
    if not isinstance(record, dict) or record.get("schema") != SCHEMA:
        raise ValueError("Unsupported RNA feature support schema")
    payload = {key: value for key, value in record.items() if key != "measurement_id"}
    if record.get("measurement_id") != content_identifier("rna_feature", payload):
        raise ValueError("RNA feature support identity mismatch")
    required = {"schema", "feature", "evidence_scope", "sampling_profile", "normalization_model",
                "geometry_scope", "values", "unavailable", "read_ids", "fragment_ids",
                "joint_only_fragment_ids", "provenance", "measurement_id"}
    if set(record) != required:
        raise ValueError("RNA support record fields differ from the schema")
    scope = record["evidence_scope"]
    if len(scope) != 2 or any(not isinstance(v, str) or not v for v in scope):
        raise ValueError("RNA evidence scope requires sample and source identities")
    feature = RNAFeature(**record["feature"])
    if record["normalization_model"] != GEOMETRY_MODEL or record["geometry_scope"] != (
            "complete_path" if feature.path_length is not None else "unbounded_local_path"):
        raise ValueError("RNA support geometry model or extent mismatch")
    profile = RNASamplingProfile.from_dict(record["sampling_profile"])
    lengths = effective_rna_feature_lengths(feature, profile)
    values = record["values"]
    expected_names = {name for unit in ("reads", "fragments") for name in (
        unit, "library_" + unit, "effective_length_" + unit, unit + "_per_million",
        unit + "_per_effective_kb", unit + "_per_effective_kb_per_million")}
    if set(values) != expected_names:
        raise ValueError("RNA support quantity fields differ from the schema")
    unavailable = {}
    for unit in ("reads", "fragments"):
        _integer(values[unit], unit)
        identifiers = record[unit[:-1] + "_ids"]
        if any(not isinstance(v, str) or not v for v in identifiers):
            raise ValueError("RNA unit identities must be nonempty strings")
        if len(identifiers) != len(set(identifiers)) or values[unit] != len(identifiers):
            raise ValueError("RNA support counts disagree with distinct unit identities")
        if values["effective_length_" + unit] != lengths[unit]:
            raise ValueError("RNA effective length disagrees with sampling geometry")
        count, denominator, length = values[unit], getattr(profile, "library_" + unit), lengths[unit]
        if values["library_" + unit] != denominator or (denominator is not None and count > denominator):
            raise ValueError("RNA library denominator mismatch")
        if length == 0 and count:
            raise ValueError("RNA support has zero modeled opportunity")
        for suffix, scale, divisors in (("per_million", 1e6, (denominator,)),
                                        ("per_effective_kb", 1e3, (length,)),
                                        ("per_effective_kb_per_million", 1e9, (denominator, length))):
            expected = None if any(v is None or v == 0 for v in divisors) else scale * count / math.prod(divisors)
            if expected is None:
                unavailable[unit + "_" + suffix] = (
                    "whole_library_denominator_unavailable" if denominator is None and "million" in suffix
                    else "sampling_geometry_unavailable" if length is None else "feature_unobservable")
            if values[unit + "_" + suffix] != expected:
                raise ValueError("RNA normalized value disagrees with counts and denominators")
    if record["unavailable"] != unavailable:
        raise ValueError("RNA unavailable reasons disagree with measurement inputs")
    joint = record["joint_only_fragment_ids"]
    if len(joint) != len(set(joint)):
        raise ValueError("Duplicate joint-only RNA template identity")
    single = values["fragments"] - len(joint)
    if not single <= values["reads"] <= 2 * single:
        raise ValueError("RNA read and template support cardinalities disagree")
    if profile.geometries and all(len(g.segments) == 1 for g in profile.geometries):
        if values["reads"] != values["fragments"] or joint:
            raise ValueError("Single-end read and fragment support must agree")
    if not set(record["joint_only_fragment_ids"]) <= set(record["fragment_ids"]):
        raise ValueError("Joint-only RNA templates are outside the support set")
    return record
