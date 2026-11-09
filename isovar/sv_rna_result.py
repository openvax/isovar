"""Assemble SV result evidence from shared observed RNA paths."""

from __future__ import annotations

from collections import Counter
from dataclasses import asdict
from hashlib import sha256

from .dna import reverse_complement_dna
from .fusion import FusionBlock
from .rna_path_geometry import PathSpans, path_breaks, position_runs, reverse_position
from .sv_rna_event import is_forward_splice
from .sv_rna_relations import LINKAGE_STATUSES


def _junction_summary(direct, segments, kinds, left, right, store):
    """Candidate-independent junction evidence, shared after discovery is frozen."""
    evidence = {o.key: o for o, *_ in direct}
    unbuilt = segments - {o.identity for o in evidence.values()}
    junction_sequences = Counter(o.sequence[a + 1:b] for o, a, b, _ in direct)
    if unbuilt:
        kinds |= store.gap_kinds.get((left, right), set()) | store.gap_kinds.get((reverse_position(right), reverse_position(left)), set())
        junction_sequences.update(inserted for identity, inserted in store.gaps.get((left, right), {}).items()
                                  if identity in unbuilt)
        junction_sequences.update(reverse_complement_dna(inserted) for identity, inserted in
                                  store.gaps.get((reverse_position(right), reverse_position(left)), {}).items() if identity in unbuilt)
    equivalent = {}
    for o in evidence.values():
        equivalent.setdefault(o.comparison_key, o)
    return dict(evidence=evidence, spanning=list(equivalent.values()), kinds=sorted(kinds),
                direct_support=store.cell_umi.support(segments),
                direct_read_lineage=store.lineage.support(segments),
                direct_observations=sorted(evidence),
                direct_junction_sequences=[list(item) for item in sorted(
                    junction_sequences.items(), key=lambda item: (-item[1], item[0]))])


def build_sv_path_result(sequence, positions, voters, store, event, annotated, adjacency, clip_support, frame_evidence,
                 orf_evidence, junction_cache=None, linked_intervals=None):
    """Serialize one path with junction, sequence and frame evidence.

    Parameters
    ----------
    sequence : str
        Observed RNA bases in path orientation.
    positions : sequence of Position or None
        One placement per observed query base; unknown placements remain None.
    voters : iterable[str]
        Original observation keys contributing to this discovery path.
    store : SvObservationStore
        Original observations, labels, lineage and join lookup.
    event : OrientedAdjacency
        Nominated adjacency policy.
    annotated : set of tuple
        Known reference joins.
    adjacency : sequence
        Direct original observations of the nominated junction.
    clip_support : mapping
        Original clipped observations by retained partner.
    frame_evidence, orf_evidence : callable
        Frame and ORF adapters consuming observed path evidence.
    junction_cache : dict, optional
        Candidate-independent summaries for one frozen observation scope.
    linked_intervals : dict, optional
        Streamed half-open co-observed path intervals for those summaries.

    Returns
    -------
    tuple[dict, dict]
        Native path result and the original observations it cites.

    Notes
    -----
    Direct support for a junction is every segment whose own observed path
    makes that join, whatever its other bases (so noisy long reads count).
    All assignments of the nominated adjacency's junction bases support its
    breakpoint junction. Counts include single records not yet built, whose
    CIGAR makes the join; built observations are cited individually.
    Returns the result and the observations it cites.
    """
    junction_cache = {} if junction_cache is None else junction_cache
    span = PathSpans(sequence, positions)
    rows = []
    for i, j in path_breaks(positions):
        left, right = positions[i], positions[j]
        relation, assignment = event.relation(left, right, j - i - 1)
        is_annotated = j == i + 1 and (left, right) in annotated
        if relation is None and is_annotated:
            continue  # Ordinary annotated splicing.
        key = (left, right, relation)
        if key not in junction_cache:
            if relation == "breakpoint_junction":
                direct, segments = adjacency, {o.identity for o, *_ in adjacency}
            else:
                direct, segments = store.joins.get((left, right), []), store.join_segments(left, right)
            kinds = {kind for *_, kind in direct}
            if segments - {o.identity for o, *_ in direct}:
                kinds |= store.gap_kinds.get((left, right), set()) | store.gap_kinds.get((reverse_position(right), reverse_position(left)), set())
            junction_cache[key] = _junction_summary(direct, segments, kinds, left, right, store)
        summary = junction_cache[key]
        # Apply the same policy to bounded discovery and streamed support.
        # CIGAR deletions are reported only as placements of the adjacency.
        # An assembled join without direct witnesses has unknown kind; absence
        # of a crossing read does not make the join a CIGAR deletion.
        kinds = set(summary['kinds'])
        if (not kinds and relation is None) or (
                kinds == {"D"} and not event.relation(left, right, j - i - 1, kind="D")[0]):
            continue
        rows.append((i, j, relation or "regional_novel_junction", assignment, is_annotated, summary))
    for i, j, assignment in event.clips(positions, len(sequence)):
        side = "donor" if assignment[0] == 0 else "acceptor"
        key = ('clip', side)
        if key not in junction_cache:
            direct = clip_support[side]
            junction_cache[key] = _junction_summary(
                direct, {o.identity for o, *_ in direct}, {"clip"}, None, None, store)
        rows.append((i, j, "breakpoint_clip_partner_unplaced", assignment, False, junction_cache[key]))
    junctions, evidence = [], {}
    for i, j, relation, assignment, is_annotated, summary in sorted(rows, key=lambda row: row[:2]):
        evidence.update(summary['evidence'])
        # Bases co-observed with the join in single reads; annotated splices are context only.
        spans = ([] if is_annotated else linked_intervals.get(id(summary), []) if linked_intervals is not None
                 else [q for q in (span(o) for o in summary['spanning']) if q])
        clip = relation == "breakpoint_clip_partner_unplaced"
        junctions.append(dict(
            # Half-open: the unplaced bases between the partners (empty for a direct join).
            query_interval=[i + 1, j], left=positions[i] and list(positions[i]),
            right=positions[j] and list(positions[j]), unplaced_bases=sequence[i + 1:j],
            forward_splice_geometry=bool(positions[i] and positions[j] and is_forward_splice(positions[i], positions[j])),
            kinds=summary['kinds'], annotated=is_annotated, relation=relation, breakpoint_assignment=assignment,
            direct_support=summary['direct_support'],
            direct_read_lineage=summary['direct_read_lineage'],
            direct_observations=summary['direct_observations'],
            direct_junction_sequences=None if clip else summary['direct_junction_sequences'],
            linked_interval=[min(min(q) for q in spans), max(max(q) for q in spans)] if spans else None))
    departures = {q: j["relation"] for j in junctions if not j["annotated"]
                  for q in range(j["query_interval"][0],
                                 len(sequence) if j["right"] is None else j["query_interval"][1] + 1)}
    frame = frame_evidence(sequence, positions, departures)
    relations = {j["relation"] for j in junctions if not j["annotated"]}
    status = next((s for s in LINKAGE_STATUSES if s in relations), "no_novel_junction")
    voting = [store.observations[key] for key in sorted(voters)]
    evidence.update((o.key, o) for o in voting)
    segments = {o.identity for o in voting}
    return dict(
        path_id=sha256(repr((sequence, positions)).encode()).hexdigest()[:24],
        sequence=sequence,
        blocks=[asdict(FusionBlock(a, b, contig, low, high, strand)) for a, b, contig, low, high, strand in position_runs(positions)],
        unplaced_intervals=_intervals(q for q, p in enumerate(positions) if p is None),
        junctions=junctions,
        event_linkage=dict(status=status, junctions=[k for k, j in enumerate(junctions)
                                                     if j["relation"] == status and not j["annotated"]],
                           somatic_causation_proven=False),
        sequence_evidence=dict(
            voting_support=store.cell_umi.support(segments),
            voting_observations=sorted(voters),
            missing_quality_reads=len({o.identity for o in voting if o.missing_qualities}),
            secondary_reads=len({o.identity for o in voting if o.secondary}),
            assembly_phase="hypothesis_not_proven_long_range_phase"),
        exploratory_orfs=orf_evidence(sequence, positions, junctions, evidence),
        **frame), evidence


def _intervals(indices):
    result = []
    for q in indices:
        if result and result[-1][1] == q:
            result[-1][1] = q + 1
        else:
            result.append([q, q + 1])
    return result
