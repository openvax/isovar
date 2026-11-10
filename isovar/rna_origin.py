"""Conservative, panel-relative nuclear/mitochondrial RNA origin evidence."""

from collections import Counter, defaultdict
from dataclasses import asdict, dataclass
import math
from tempfile import TemporaryDirectory

import edlib

from .dna import reverse_complement_dna
from .read_collector import ReadCollector
from .read_identity import segment_identity
from .rna_evidence import content_identifier, evidence_set, rna_support
from .rna_feature_measurement import rna_path_offsets
from .rna_feature_support import RNAFeature
from .rna_library import original_library_context
from .rna_observation import ObservationGeometryCache, build_segment_observations
from .rna_path_geometry import PathPositions, position_runs
from .rna_reconstruction import placement_uncertainty


SCHEMA = "isovar.rna_origin.v1"
ORIGINS = frozenset(("nuclear", "mitochondrial"))
STATUSES = ("nuclear_favored", "mitochondrial_favored", "ambiguous", "conflicting", "unassessed")
_MASK_EQUALITIES = [("N", base) for base in "ACGT"]


@dataclass(frozen=True)
class RNAOriginPolicy:
    """Uncalibrated context-comparison thresholds, recorded in every result.

    Q20 and 30 reliable flank bases are required by default. The best origin
    must fit within 5% context edits and beat the other origin by three edits,
    plus one extra edit for each masked context base. Missing QUAL cannot
    establish an origin. No MAPQ, allele fraction or protein score is a vote.
    """

    min_base_quality: int = 20
    min_context_bases: int = 30
    max_edit_fraction: float = 0.05
    min_edit_margin: int = 3

    def __post_init__(self):
        for name, minimum in (("min_base_quality", 0), ("min_context_bases", 1), ("min_edit_margin", 1)):
            value = getattr(self, name)
            if type(value) is not int or value < minimum:
                raise ValueError(f"{name} must be an integer >= {minimum}")
        if (isinstance(self.max_edit_fraction, bool) or not isinstance(self.max_edit_fraction, (int, float))
                or not math.isfinite(self.max_edit_fraction) or not 0 <= self.max_edit_fraction < 1):
            raise ValueError("max_edit_fraction must be finite and in [0, 1)")


@dataclass(frozen=True)
class RNAOriginHypothesis:
    """An explicitly labeled nucleotide path carrying the nominated target.

    ``positions`` uses the shared (contig, zero-based base, strand) RNA-path
    convention; inserted bases are None. ``target_interval`` is half-open,
    nonempty and includes placed anchors for an insertion/deletion junction.
    All competitors must carry the same target sequence up to reverse
    complement. Supply both origin groups and sufficient flanking context;
    labels and panel completeness are caller declarations, not discoveries.
    """

    identifier: str
    origin: str
    sequence: str
    positions: tuple
    target_interval: tuple[int, int]

    def __post_init__(self):
        if not isinstance(self.identifier, str) or not self.identifier or self.origin not in ORIGINS:
            raise ValueError("Origin hypothesis requires an identifier and nuclear/mitochondrial origin")
        if not isinstance(self.sequence, str) or not self.sequence or set(self.sequence.upper()) - set("ACGT"):
            raise ValueError("Origin hypothesis requires resolved ACGT sequence")
        object.__setattr__(self, "sequence", self.sequence.upper())
        positions = PathPositions(None if p is None else tuple(p) for p in self.positions)
        if len(positions) != len(self.sequence) or any(p is not None and (
                len(p) != 3 or not isinstance(p[0], str) or not p[0] or type(p[1]) is not int
                or p[1] < 0 or p[2] not in ("+", "-")) for p in positions):
            raise ValueError("Origin sequence requires matching per-base RNA path positions")
        interval = tuple(self.target_interval)
        RNAFeature(self.identifier, self.identifier, (interval,), "compound", path_length=len(self.sequence))
        if not any(p is not None for p in positions[interval[0]:interval[1]]):
            raise ValueError("Origin target requires placed anchor bases")
        object.__setattr__(self, "positions", positions)
        object.__setattr__(self, "target_interval", interval)

    @property
    def target_sequence(self):
        return self.sequence[slice(*self.target_interval)]

    def to_dict(self):
        """Return a JSON-compatible path declaration, without inferred origin."""
        return dict(identifier=self.identifier, origin=self.origin, sequence=self.sequence,
                    positions=list(self.positions), target_interval=list(self.target_interval))


def _target_views(observations, hypotheses):
    """Locate the target only on actual, coordinate-compatible observations."""
    views = {}
    for observation in observations:
        for hypothesis in hypotheses:
            a, b = hypothesis.target_interval
            for offset in rna_path_offsets(hypothesis.positions, observation):
                lo, hi = a - offset, b - offset
                if (0 <= lo < hi <= len(observation.sequence)
                        and observation.positions[lo:hi] == hypothesis.positions[a:b]):
                    # Complementary targets are allowed for opposite-oriented
                    # models. At the same placement they remain alternatives,
                    # rather than letting panel order erase a matching allele.
                    key = (observation.key, lo, hi, hypothesis.target_sequence)
                    views[key] = (observation, lo, hi,
                                  observation.sequence[lo:hi] == hypothesis.target_sequence)
    return list(views.values())


def _orientations(hypothesis, target):
    a, b = hypothesis.target_interval
    sequence = hypothesis.sequence
    if hypothesis.target_sequence == target:
        yield sequence[:a], sequence[b:]
    if reverse_complement_dna(hypothesis.target_sequence) == target:
        yield reverse_complement_dna(sequence[b:]), reverse_complement_dna(sequence[:a])


def _distance(query, target):
    if not query:
        return 0
    return edlib.align(query, target, mode="SHW", task="distance",
                       additionalEqualities=_MASK_EQUALITIES)["editDistance"]


def _compare_context(observation, lo, hi, hypotheses, policy):
    result = dict(status="unassessed", reason=None, comparisons=[])
    sequence, qualities = observation.sequence, observation.qualities
    if len(qualities) != len(sequence) or any(q is None for q in qualities):
        return dict(result, reason="missing_base_qualities")
    if any(q < policy.min_base_quality for q in qualities[lo:hi]):
        return dict(result, reason="low_quality_target")
    masked = "".join(base if base in "ACGT" and q >= policy.min_base_quality else "N"
                     for base, q in zip(sequence, qualities))
    left, right = masked[:lo], masked[hi:]
    unreliable = (left + right).count("N")
    reliable = len(left) + len(right) - unreliable
    result.update(reliable_context_bases=reliable, masked_context_bases=unreliable)
    if reliable < policy.min_context_bases:
        return dict(result, reason="insufficient_reliable_context")
    incomplete = False
    for hypothesis in hypotheses:
        costs = []
        for before, after in _orientations(hypothesis, sequence[lo:hi]):
            if len(before) < len(left) or len(after) < len(right):
                incomplete = True
                continue
            costs.append(_distance(left[::-1], before[::-1]) + _distance(right, after))
        result["comparisons"].append(dict(hypothesis_id=hypothesis.identifier, origin=hypothesis.origin,
                                          edit_distance=min(costs) if costs else None))
    if incomplete:
        return dict(result, reason="insufficient_panel_context")
    scores = {origin: min(c["edit_distance"] for c in result["comparisons"] if c["origin"] == origin)
              for origin in sorted(ORIGINS)}
    best = min(scores.values())
    margin = abs(scores["nuclear"] - scores["mitochondrial"])
    result.update(origin_edit_distances=scores, edit_margin=margin,
                  required_edit_margin=policy.min_edit_margin + unreliable)
    if best > policy.max_edit_fraction * reliable:
        return dict(result, reason="poor_fit_to_panel")
    if margin < policy.min_edit_margin + unreliable:
        return dict(result, status="ambiguous", reason="origins_not_separated")
    return dict(result, status=min(scores, key=scores.get) + "_favored", reason="context_edit_margin")


def _combine(statuses):
    """Alternatives and mates are constraints, never independent majority votes."""
    statuses = set(statuses)
    favored = statuses & {"nuclear_favored", "mitochondrial_favored"}
    if "conflicting" in statuses or len(favored) > 1:
        return "conflicting"
    if len(statuses) == 1:
        return next(iter(statuses))
    return "ambiguous"


def _sample_assessment(states):
    """Independent fragments: uninformative ones are counted but cannot veto.

    Opposing favored fragments stay mixed and any conflicting fragment stays
    conflicting. Otherwise one favored origin summarizes the informative
    fragments; ambiguous/unassessed fragments remain in ``support_by_status``.
    """
    favored = states & {"nuclear_favored", "mitochondrial_favored"}
    if not states:
        return "no_target_support"
    if len(favored) > 1:
        return "mixed_origins"
    if "conflicting" in states:
        return "conflicting"
    if favored:
        return next(iter(favored))
    return "ambiguous" if "ambiguous" in states else "unassessed"


def assess_rna_origin(bam, hypotheses, *, sample_id, source, policy=None,
                      read_collector=None, library_index=None, reference_name=None,
                      read_filter_identity=None, scratch_dir=None):
    """Assess nuclear versus mitochondrial context of target-bearing RNA.

    Parameters
    ----------
    bam : pysam.AlignmentFile
        Original library, or explicitly limited input. Its mapped target-bearing
        segments are selected through the shared original-library index.
    hypotheses : iterable of RNAOriginHypothesis
        Explicit competitors carrying the same allele/anchored target. At least
        one nuclear and one mitochondrial hypothesis are required.
    sample_id, source : str
        Nonempty sample and original-library identities; libraries never pool.
    policy : RNAOriginPolicy, optional
        Auditable, uncalibrated comparison thresholds.
    read_collector : ReadCollector, optional
        Shared RNA-path eligibility policy (retains MAPQ-zero and secondary
        placements by default). Missing QUAL remains visible but unassessed.
    library_index : RnaLibraryIndex, optional
        Reuse a validated original-library index without rescanning the BAM.
    reference_name, read_filter_identity : str, optional
        Reference/filter identity for original-library index validation.
    scratch_dir : path, optional
        Parent of the disposable index directory when no index is supplied.

    Returns
    -------
    dict
        Native JSON-compatible result with per-fragment/placement assessments,
        disjoint fragment categories, scoped supporting read sets, policy and
        input/panel provenance. Favored origins are relative to the supplied
        panel, not proofs of origin or calibrated probabilities. Translation
        tables, candidate admission and RNA counts elsewhere are unchanged.
    """
    hypotheses = tuple(hypotheses)
    if not hypotheses or any(not isinstance(h, RNAOriginHypothesis) for h in hypotheses):
        raise ValueError("Origin assessment requires RNAOriginHypothesis records")
    hypotheses = tuple(sorted(hypotheses, key=lambda h: h.identifier))
    if len({h.identifier for h in hypotheses}) != len(hypotheses):
        raise ValueError("Origin hypothesis identifiers must be unique")
    if {h.origin for h in hypotheses} != ORIGINS:
        raise ValueError("Origin panel requires nuclear and mitochondrial hypotheses")
    targets = {min(h.target_sequence, reverse_complement_dna(h.target_sequence)) for h in hypotheses}
    if len(targets) != 1:
        raise ValueError("Origin hypotheses must carry the same target allele/sequence")
    if any(not isinstance(s, str) or not s for s in (sample_id, source)):
        raise ValueError("Origin evidence requires nonempty sample_id and source")
    policy = RNAOriginPolicy() if policy is None else policy
    if not isinstance(policy, RNAOriginPolicy):
        raise ValueError("policy must be an RNAOriginPolicy")
    lengths = dict(zip(bam.references, bam.lengths))
    regions = []
    for h in hypotheses:
        for _, _, contig, a, b, _ in position_runs(h.positions):
            if contig not in lengths or b > lengths[contig]:
                raise ValueError("Origin path is outside the supplied alignment reference")
            regions.append((contig, a, b))
    scope = [sample_id, source]
    collector = read_collector or ReadCollector.for_rna_paths()
    cache, rejections, uncertainty = ObservationGeometryCache(), Counter(), Counter()
    categories, evidence_sets, fragments = {s: set() for s in STATUSES}, {}, []

    def support(identities):
        record = evidence_set(scope, identities)
        evidence_sets[record["evidence_set_id"]] = record
        return dict(rna_support(identities), evidence_set_id=record["evidence_set_id"])

    with TemporaryDirectory(prefix="isovar-origin-", dir=scratch_dir) as directory:
        with original_library_context(bam, collector, directory, library_index=library_index,
                reference_name=reference_name, read_filter_identity=read_filter_identity) as groups:
            for template in groups.templates(regions):
                uncertainty.update(placement_uncertainty(template))
                records = defaultdict(list)
                for record in template:
                    records[segment_identity(record)].append(record)
                segments, keys = [], set()
                for identity, members in sorted(records.items()):
                    observations = build_segment_observations(identity, members, collector, rejections, {}, cache)
                    views = _target_views(observations, hypotheses)
                    if not any(present for _, _, _, present in views):
                        continue
                    assessments = []
                    for observation, lo, hi, present in views:
                        assessed = (_compare_context(observation, lo, hi, hypotheses, policy) if present else
                                    dict(status="ambiguous", reason="conflicting_target_placement", comparisons=[]))
                        assessments.append(dict(assessed, observation_id=content_identifier(
                            "origin_observation", [scope, identity, observation.key, lo, hi]),
                            original_record_ids=list(observation.records)))
                    status = _combine(a["status"] for a in assessments)
                    keys.add(identity)
                    segments.append(dict(read_id=content_identifier("segment", [scope, identity]),
                                         status=status, observations=assessments))
                if segments:
                    status = _combine(s["status"] for s in segments)
                    categories[status].update(keys)
                    fragments.append(dict(fragment_id=content_identifier("fragment", [scope, next(iter(keys))[:2]]),
                                          status=status, rna_support=support(keys), segments=segments))
            provenance = dict(input_identity=groups.prepared.identity, input_extent="declared_input",
                              read_policy=groups.prepared.identity["eligibility"],
                              observation_rejections=dict(rejections), placement_uncertainty=dict(uncertainty))
    status = _sample_assessment({f["status"] for f in fragments})
    panel = [h.to_dict() for h in hypotheses]
    return dict(schema=SCHEMA, assessment=status, evidence_scope=scope, policy=asdict(policy),
                comparison_policy="allele_anchored_origin_context.v1", hypotheses=panel,
                panel_id=content_identifier("origin_panel", panel), fragments=fragments,
                rna_support=support(set().union(*categories.values())),
                support_by_status={s: support(keys) for s, keys in categories.items()},
                evidence_sets=evidence_sets, provenance=provenance,
                limitations=["Panel-relative heuristic, not calibrated origin probability or NuMT exclusion.",
                             "Only actual target-bearing segments establish allele-to-context linkage; disjoint mates do not vote.",
                             "Unmapped reads and missing nuclear/patient-specific paths are not discovered by this assessment.",
                             "Library chimeras can mimic linkage; no independent DNA validation is implied.",
                             "No translation-table change, candidate filtering, DNA heteroplasmy or abundance inference."])
