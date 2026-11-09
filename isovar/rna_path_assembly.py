"""Exact RNA overlap extension with original-segment placement constraints."""

from __future__ import annotations

from collections import defaultdict

from .dna import reverse_complement_dna
from .rna_observation import (
    ObservationChoices,
    RnaObservationStore,
    group_observation_choices,
)
from .rna_path_geometry import exact_overlap_agrees, reverse_path, reverse_position

_LOCAL = 20  # Placements this close are local (base/indel) alternatives, not splice choices.


SEQUENCE_BRIDGE_BASES = 12  # Exact clipped seed, followed by the usual full overlap check.


def _sequence_overlapping(index, sequence, positions, need):
    """Find bounded exact overlaps through an unplaced terminal RNA sequence.

    A 12-base clipped suffix is only a lookup seed. The complete shared span
    must meet ``need`` and agree base by base, including every jointly placed
    coordinate. A sequence-only bridge does not assign its clipped bases to a
    genomic locus or establish the nominated DNA event's cause.
    """
    n = len(sequence)
    if n < SEQUENCE_BRIDGE_BASES or any(p is not None for p in positions[-SEQUENCE_BRIDGE_BASES:]):
        return []
    motif = sequence[-SEQUENCE_BRIDGE_BASES:]
    if isinstance(index, RnaObservationStore):
        for identity in index.groups:
            index.build(identity)
    matches = []
    for o in index.observations.values():
        start = o.sequence.find(motif)
        while start >= 0:
            offset = n - SEQUENCE_BRIDGE_BASES - start
            if (offset <= n - need and offset + len(o.sequence) > n
                    and exact_overlap_agrees(sequence, positions, o, offset, max(0, offset), n)):
                matches.append((o, offset))
            start = o.sequence.find(motif, start + 1)
    groups = group_observation_choices(matches)
    cap = index.cap if isinstance(index, RnaObservationStore) else len(index.observations)
    if len(groups) > cap:
        if isinstance(index, RnaObservationStore):
            index.limited = True
        groups.sort(key=lambda g: (-(n - max(0, g[1])), -len(g[2]),
                                   -(g[1] + len(g[0].sequence) - n), g[0].key))
        groups = groups[:cap]
    return groups


def supported_alternatives(support, min_fragments, min_fraction, min_local_fraction):
    """Apply explicit fragment thresholds to competing continuation bases.

    Parameters
    ----------
    support : mapping
        Nonempty (base, placement) to original-fragment count mapping.
    min_fragments : int
        Minimum original fragments for a retained alternative.
    min_fraction : float
        Minimum count relative to the best-supported alternative.
    min_local_fraction : float
        Additional relative threshold for unplaced or nearby placements
        (within 20 bases on the same contig and strand).

    Returns
    -------
    set
        Retained (base, placement) alternatives; empty if none meets the policy.
    """
    best_key = max(support, key=lambda key: (support[key], repr(key)))
    best = support[best_key]

    def local(key):
        a, b = key[1], best_key[1]
        return a is None or b is None or (a[::2] == b[::2] and abs(a[1] - b[1]) <= _LOCAL)

    return {key for key, count in support.items()
            if count >= min_fragments and count >= min_fraction * best
            and (key == best_key or not local(key) or count >= min_local_fraction * best)}


def _vote_members(groups):
    return ((observation, offset) for _, offset, members in groups for observation in members.values())


def _gather_extenders(sequence, positions, index, need, used, carried, conflicts=None):
    """Compare each sequence/placement once, including already agreeing groups."""
    n = len(sequence)
    found = {(o.comparison_key, offset): (o, offset, members) for o, offset, members in carried
             if offset <= n - need and offset + len(o.sequence) > n}
    grouped = index.grouped_overlapping(positions, n - need, n)
    if getattr(index, 'sequence_overlaps', False):
        grouped += _sequence_overlapping(index, sequence, positions, need)
    for o, offset, members in grouped:
        if offset > n - need or offset + len(o.sequence) <= n:
            continue
        key = o.comparison_key, offset
        previous = found.get(key)
        keys = members.keys() - previous[2].keys() if previous else members.keys()
        if not keys or (previous is None and not exact_overlap_agrees(sequence, positions, o, offset, max(0, offset), n)):
            continue
        permitted = ObservationChoices((k, members[k]) for k in keys
                                   if used.get(members[k].identity, k) == k)
        if conflicts is not None:
            conflicts.update(members[key].identity for key in keys - permitted.keys())
        if permitted:
            if previous is not None:
                permitted = ObservationChoices(previous[2] | permitted)
            found[key] = o, offset, permitted
    return list(found.values())


def _refine_unplaced(positions, extenders):
    """Assign a clipped base only when all recruited reads agree on its locus.

    The sequence bridge is still an assembly hypothesis. Refining coordinates
    lets a reconstructed ordinary splice compete with an SV interpretation;
    conflicting placements stay unplaced rather than choosing one arbitrarily.
    """
    updated = list(positions)
    occupied = {p for p in positions if p is not None}
    for q, current in enumerate(positions):
        if current is not None:
            continue
        choices = {o.positions[q - offset] for o, offset, _ in extenders
                   if offset <= q < offset + len(o.positions) and o.positions[q - offset] is not None}
        if len(choices) == 1:
            position, = choices
            if position not in occupied:
                updated[q] = position
                occupied.add(position)
    return tuple(updated)


def _choose_voters(used, groups, rejected):
    """Preserve the original smallest-key tie-break between alternate observations."""
    choices = {}
    for _, _, members in groups:
        if rejected and not members.keys().isdisjoint(rejected):
            members = ObservationChoices((k, o) for k, o in members.items() if k not in rejected)
        incoming = members.choices
        common = choices.keys() & incoming.keys()
        first = {identity: min(choices[identity], incoming[identity]) for identity in common}
        choices.update(incoming)
        choices.update(first)
    choices.update(used)  # A previously selected placement cannot switch.
    return choices


def _carry_votes(groups, used):
    allowed = set(used.values())
    carried = []
    for o, offset, members in groups:
        keys = members.keys() & allowed
        keys.update(k for k in members.keys() - keys if members[k].identity not in used)
        if keys:
            carried.append((o, offset, members if len(keys) == len(members) else
                            ObservationChoices((k, members[k]) for k in keys)))
    return carried


def _votes(active, t):
    """Group reads by their (base, placement) at path offset ``t``.

    An unplaced base agrees with every placement of the same base, so a
    soft-clipped or inserted copy of a base does not split its placed support.
    """
    votes = defaultdict(list)
    for observation, offset, members in active:
        if t - offset < len(observation.sequence):
            votes[observation.sequence[t - offset], observation.positions[t - offset]].append(
                (observation, offset, members))
    for base, position in [key for key in votes if key[1] is None]:
        placed = [key for key in votes if key[0] == base and key[1] is not None]
        if placed:
            for key in placed:
                votes[key].extend(votes[base, None])
            del votes[base, None]
    return votes


def _extend(sequence, positions, index, min_overlap, thresholds, budget, pruned, notes, mirrored,
            initial_used=None, conflicts=None):
    """Extend 3' through exact overlaps which share a placed base with the path.

    Continuations are added while every agreeing read agrees on the next base
    and its placement. At a disagreement, weakly supported alternatives are
    pruned (and reported) only when another alternative meets both
    thresholds; otherwise every alternative becomes its own path. Support is
    counted over every read overlapping the path end. A branch keeps the reads
    which chose it, so it can continue through unplaced sequence (e.g.
    soft-clipped tails) where no new read can be anchored. Each finished
    path is returned with the keys of the observations which voted for it and
    why this end stopped: ``observations_end`` when no observation extends it,
    ``repeated_genomic_position`` or ``path_limit``.
    """
    done, stack = [], [(sequence, positions, dict(initial_used or {}), ())]
    while stack:
        sequence, positions, used, carried = stack.pop()
        n = len(sequence)
        need = min(min_overlap, n)
        extenders = _gather_extenders(sequence, positions, index, need, used, carried, conflicts)
        if not extenders:
            done.append((sequence, positions, frozenset(used.values()), "observations_end"))
            continue
        if getattr(index, 'sequence_overlaps', False):
            positions = _refine_unplaced(positions, extenders)
        placed = {p for p in positions if p is not None}
        bases, added, active, forks, stopped, dropped = [], [], extenders, None, False, set()
        extender_count = sum(len(members) for _, _, members in extenders)
        while True:
            votes = _votes(active, n + len(bases))
            if not votes or (bases and (len(votes) > 1 or
                                       2 * sum(len(g[2]) for v in votes.values() for g in v) < extender_count)):
                # Reads ending within this step leave fewer voters: decide a
                # disagreement, or continue once most extenders have ended,
                # only after re-collecting every read overlapping the new end
                # (so one long read cannot decide the rest of the path alone).
                break
            if len(votes) > 1:
                support = {key: len(set().union(*(g[2].fragments for g in v))) for key, v in votes.items()}
                strong = supported_alternatives(support, *thresholds)
                if strong and len(strong) < len(votes):
                    anchor = next((p for p in reversed(positions + tuple(added)) if p is not None), None)
                    for key in sorted(set(votes) - strong, key=repr):
                        base, position = key
                        pruned.append(dict(
                            direction="5prime" if mirrored else "3prime",
                            after=None if anchor is None else list(reverse_position(anchor) if mirrored else anchor),
                            base=reverse_complement_dna(base) if mirrored else base,
                            placement=None if position is None else list(reverse_position(position) if mirrored else position),
                            fragments=support[key], observations=sorted({o.key for o, _ in _vote_members(votes[key])})))
                        dropped.update(o.key for o, _ in _vote_members(votes[key]))
                    votes = {key: votes[key] for key in strong}
                if len(votes) > 1:
                    forks = votes
                    break
            ((base, position), active), = votes.items()
            if position is not None:
                if position in placed:
                    # A path may not revisit a placement (e.g. around a tandem
                    # duplication); copy number cannot be read from overlaps.
                    notes.add("repeated_genomic_position")
                    stopped = True
                    break
                placed.add(position)
            bases.append(base)
            added.append(position)
        sequence, positions = sequence + "".join(bases), positions + tuple(added)
        if stopped or (forks and len(done) + len(stack) + len(forks) > budget):
            if forks:
                notes.add("path_limit")
            voters = _choose_voters(used, extenders, dropped)
            done.append((sequence, positions, frozenset(voters.values()),
                         "repeated_genomic_position" if stopped else "path_limit"))
        elif forks:
            for key in sorted(forks, key=repr):
                chosen = set().union(*(g[2].keys() for g in forks[key]))
                rejected = set().union(*(g[2].keys() for other, v in forks.items() if other != key for g in v)) - chosen
                branch = _choose_voters(used, extenders, rejected | dropped)
                stack.append((sequence + key[0], positions + (key[1],), branch, _carry_votes(forks[key], branch)))
        else:
            used = _choose_voters(used, extenders, dropped)
            stack.append((sequence, positions, used, _carry_votes(active, used) if bases else ()))
    return done


def _mirrored_choices(voters, index):
    """Keep each chosen original observation on its reciprocal orientation."""
    return {index.observations[key].identity: key[:-1] + ('-' if key.endswith('+') else '+')
            for key in voters}


def extend_path_both_directions(sequence, positions, index, min_overlap, thresholds, budget, pruned, notes,
                 *, seed_observation=None):
    """Extend an observed RNA seed with consistent original-segment choices.

    Parameters
    ----------
    sequence : str
        Observed RNA bases in the orientation of positions; no reference filling.
    positions : sequence of Position or None
        One placement per query base: (contig, zero-based base, strand), or
        None for an observed base whose placement is unresolved.
    index : ObservationIndex or RnaObservationStore
        Original observations, with reciprocal views for both extension ends.
    min_overlap : int
        Minimum exact overlap, clipped to the current seed length.
    thresholds : tuple[int, float, float]
        Minimum fragments, relative support and local relative support;
        see supported_alternatives. (0, 0, 0) retains all observed choices.
    budget : int
        Maximum distinct returned paths across both extension orders.
    pruned : list
        Receives recorded alternatives rejected by the threshold policy.
    notes : set[str]
        Receives explicit limits and unresolved placement-choice notes.
    seed_observation : RnaObservation, optional
        For a whole-read seed, its actual indexed observation. Prevents an
        alternative placement of that same segment extending its own seed.

    Returns
    -------
    iterator of tuple
        (sequence, placements, voter_keys, end_states). Voters are discovery
        contributors, not full-sequence witnesses or molecule counts.

    Raises
    ------
    ValueError
        A supplied whole-read seed is not its actual indexed observation.

    Notes
    -----
    No reference sequence fills a gap. Opposite-end placement conflicts remain
    separate hypotheses, bounded by one shared path budget.
    """
    initial_choices = {}
    if seed_observation is not None:
        if (index.observations.get(seed_observation.key) is not seed_observation or
                seed_observation.comparison_key != (sequence, positions)):
            raise ValueError("A whole-path seed must be its actual indexed observation")
        initial_choices = {seed_observation.identity: seed_observation.key}
    seen, conflicts = set(), set()
    for reverse_first in (False, True):
        if reverse_first:
            if not conflicts:
                break
            notes.add('incompatible_segment_placements_retained_as_alternatives')
        seed = reverse_path(sequence, positions) if reverse_first else (sequence, positions)
        seed_choices = (_mirrored_choices(initial_choices.values(), index) if reverse_first else initial_choices)
        for first_sequence, first_positions, voters, first_end in _extend(
                *seed, index, min_overlap, thresholds, budget, pruned, notes, reverse_first,
                initial_used=seed_choices, conflicts=conflicts if seed_observation is not None else None):
            choices = _mirrored_choices(voters, index)
            for second_sequence, second_positions, more, second_end in _extend(
                    *reverse_path(first_sequence, first_positions), index, min_overlap, thresholds,
                    budget, pruned, notes, not reverse_first, initial_used=choices, conflicts=conflicts):
                path = ((second_sequence, second_positions) if reverse_first else
                        reverse_path(second_sequence, second_positions))
                if path in seen:
                    continue
                if len(seen) >= budget:
                    notes.add('path_limit')
                    return
                seen.add(path)
                ends = ({'5prime': first_end, '3prime': second_end} if reverse_first else
                        {'5prime': second_end, '3prime': first_end})
                yield (*path, voters | (more - set(choices.values())), ends)
