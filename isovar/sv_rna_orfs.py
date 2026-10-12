"""Exploratory ATG ORFs and their observed RNA evidence, not initiation calls."""

from bisect import bisect_right
from itertools import chain, islice

from .genetic_code import standard_genetic_code
from .orf_inclusion import annotate_orf_inclusion
from .orf_start import annotate_orf_start
from .sv_rna_relations import EVENT_CROSSING_RELATIONS


_PLACEMENT_INDEX = object()
_WITNESS_CACHE_ENTRIES = 4096


def _sequence_offsets(sequence, dna):
    offset = sequence.find(dna)
    while offset >= 0:
        yield offset
        offset = sequence.find(dna, offset + 1)


class _PlacementIndex:
    """Compare placed coordinates by linear runs, allowing unplaced bases."""

    def __init__(self, positions):
        # sv_rna imports the ORF enumerator; run geometry is available at call
        # time. Retaining positions also prevents recycling its cache-key id.
        from .sv_rna import _runs

        self.positions = positions
        self.runs = _runs(positions)
        self.ends = [run[1] for run in self.runs]

    def agrees(self, other, start, end, shift):
        """Every base placed on both sides has the same genomic placement."""
        i, j = bisect_right(self.ends, start), bisect_right(other.ends, start + shift)
        while i < len(self.runs) and j < len(other.runs):
            p0, p1, pc, pl, ph, ps = self.runs[i]
            q0, q1, qc, ql, qh, qs = other.runs[j]
            left, right = max(start, p0, q0 - shift), min(end, p1, q1 - shift)
            if left < right:
                if pc != qc or ps != qs:
                    return False
                own = pl + left - p0 if ps == '+' else ph - 1 - left + p0
                observed = ql + left + shift - q0 if qs == '+' else qh - 1 - left - shift + q0
                if own != observed:
                    return False
            if min(p1, q1 - shift) >= end:
                break
            if p1 <= q1 - shift:
                i += 1
            if q1 - shift <= p1:
                j += 1
        return True


def _witnesses(sequence, positions, start, end, junctions, observations, group_cache=None):
    """Exact full-interval witnesses with compatible placements at this join.

    Sequence agreement includes unplaced bases. Two partial reads never
    become a full witness. Different assignments of junction homology are
    allowed only at unplaced bases, with shared placed anchors on each
    available side. A shorter ORF may need additional sequence to link it to
    the other flank; that interval must also match in this same observation.
    """
    # One cache belongs to one frozen observation store. Quality/metadata do
    # not affect sequence matching and are read from each original afterwards.
    group_cache = {} if group_cache is None else group_cache
    dna = sequence[start:end]
    witnesses = {}
    matches, consistent_intervals = {}, {}
    cached_matches = 0

    def offsets(text):
        nonlocal cached_matches
        if text in matches:
            return iter(matches[text])
        iterator = _sequence_offsets(text, dna)
        room = _WITNESS_CACHE_ENTRIES - cached_matches
        if room <= 0:
            return iterator
        found = tuple(islice(iterator, room))
        if len(found) < room:
            matches[text] = found
            cached_matches += len(found) + 1
            return iter(found)
        # Repetitive sequences still yield every occurrence, without filling
        # an unbounded offsets table when no occurrence has valid placements.
        return chain(found, iterator)

    def index(pos):
        key = _PLACEMENT_INDEX, id(pos)
        if key not in group_cache:
            group_cache[key] = _PlacementIndex(pos)
        return group_cache[key]

    for junction in junctions:
        # Flanking placed bases around the junction's unplaced interval.
        left, right = junction["query_interval"][0] - 1, junction["query_interval"][1]
        keys = tuple(junction["direct_observations"])
        if keys not in group_cache:
            groups = {}
            for key in keys:
                observation = observations[key]
                groups.setdefault(observation.comparison_key, []).append(observation)
            group_cache[keys] = list(groups.values())
        for members in group_cache[keys]:
            observation = members[0]
            comparison = observation.comparison_key
            for offset in offsets(observation.sequence):
                shift = offset - start
                lo, hi = max(0, -shift), min(len(sequence), len(observation.sequence) - shift)

                def shared(q):
                    return positions[q] is not None and positions[q] == observation.positions[q + shift]

                anchors = []
                if junction["left"] is not None:
                    anchors.append(next((q for q in range(min(left, end - 1), lo - 1, -1) if shared(q)), None))
                if junction["right"] is not None:
                    anchors.append(next((q for q in range(max(right, start), hi) if shared(q)), None))
                if anchors and all(q is not None for q in anchors):
                    a, b = min(start, *anchors), max(end, *(q + 1 for q in anchors))
                    interval = comparison, shift, a, b
                    if interval in consistent_intervals:
                        consistent = consistent_intervals[interval]
                    else:
                        consistent = (
                            sequence[a:b] == observation.sequence[a + shift:b + shift]
                            and index(positions).agrees(index(observation.positions), a, b, shift))
                        if len(consistent_intervals) < _WITNESS_CACHE_ENTRIES:
                            consistent_intervals[interval] = consistent
                else:
                    consistent = False
                if consistent:
                    link_interval = [a + shift, b + shift]
                    for original in members:
                        a = (original.query_interval[1] - offset - len(dna) if original.reverse_complement
                             else original.query_interval[0] + offset)
                        witness = witnesses.setdefault((original.key, offset), dict(
                            observation=original.key, observation_interval=[offset, offset + len(dna)],
                            original_query_interval=[a, a + len(dna)], junction_links=[]))
                        link = dict(junction_query_interval=list(junction["query_interval"]),
                                    observation_interval=list(link_interval))
                        if link not in witness["junction_links"]:
                            witness["junction_links"].append(link)
    return [witnesses[key] for key in sorted(witnesses)]


def _boundaries(junction):
    """Named RNA boundaries; unplaced bases need not be a DNA insertion."""
    start, end = junction["query_interval"]
    if start == end and junction["left"] is not None and junction["right"] is not None:
        return [("flank_to_flank", end)]
    return ([("donor_to_unplaced", start)] if junction["left"] is not None else []) + (
        [("unplaced_to_acceptor", end)] if junction["right"] is not None else [])


def exploratory_orfs(sequence, positions, junctions, observations, references, rna_support,
                     min_amino_acids, max_candidates, lineage=None, competing_splices=(), inclusion_thresholds=None,
                     witness_groups=None):
    """Enumerate bounded ATG candidates crossing an event-related RNA join.

    These are potential translations of observed sequence, separate from
    annotated CDS frame transfer. Context is reported, never scored as proof
    of initiation. A missing stop means the retained sequence ends first.
    ``references`` are the supplied ``FusionReference`` transcript models.
    ``inclusion_thresholds`` holds keyword arguments for
    ``annotate_orf_inclusion``; None uses its defaults.
    """
    related = [(i, j, _boundaries(j)) for i, j in enumerate(junctions)
               if not j["annotated"] and j["relation"] in EVENT_CROSSING_RELATIONS]
    candidates, limited = [], False
    # Next in-frame stop for each ATG, without translating every long suffix.
    next_stop = [None, None, None]
    starts = []
    for q in range(len(sequence) - 3, -1, -1):
        codon = sequence[q:q + 3]
        if codon in standard_genetic_code.stop_codons:
            next_stop[q % 3] = q
        elif codon == "ATG":
            stop = next_stop[q % 3]
            coding_end = stop if stop is not None else q + 3 * ((len(sequence) - q) // 3)
            end = coding_end + (3 if stop is not None else 0)
            crossed = [(i, j, [(kind, b) for kind, b in boundaries if q < b < end])
                       for i, j, boundaries in related if any(q < b < end for _, b in boundaries)]
            if crossed and (coding_end - q) // 3 >= min_amino_acids:
                starts.append((q, coding_end, stop is not None, crossed))
    limited = len(starts) > max_candidates
    for start, coding_end, has_stop, crossed in reversed(starts[-max_candidates:]):
        end = coding_end + (3 if has_stop else 0)
        aa, _ = standard_genetic_code.translate(sequence[start:coding_end], first_codon_is_start=True)
        witnesses = _witnesses(sequence, positions, start, end, [j for _, j, _ in crossed], observations,
                              group_cache=witness_groups)
        segments = {observations[w["observation"]].identity for w in witnesses}
        start_evidence = annotate_orf_start(sequence, positions, start, references)
        start_evidence = annotate_orf_inclusion(
            start_evidence, sequence, positions, end, references, observations, witnesses,
            competing_splices=competing_splices, **(inclusion_thresholds or {}))
        for assessment in start_evidence["assessments"]:
            inclusion = assessment.get("splice_inclusion")
            if inclusion is not None:
                qualified = {observations[key].identity for key in inclusion["qualified_observations"]}
                inclusion["support"] = rna_support(qualified)
                inclusion["read_lineage"] = lineage(qualified) if lineage is not None else None
        comparisons = start_evidence["reference_comparisons"]
        minus_three = sequence[start - 3] if start >= 3 else None
        plus_four = sequence[start + 3] if start + 3 < len(sequence) else None
        candidates.append(dict(
            query_interval=[start, end], amino_acids=aa, ends_with_stop_codon=has_stop,
            initiation_observed=False, translation_observed=False,
            start_codon_basis="observed_ATG", annotated_start=any(c["start_kind"] == "annotated_CDS_start" for c in comparisons),
            start_context=dict(sequence=sequence[max(0, start - 6):min(len(sequence), start + 9)],
                               query_interval=[max(0, start - 6), min(len(sequence), start + 9)],
                               minus_three=minus_three, plus_four=plus_four),
            reference_comparisons=comparisons, start_evidence=start_evidence,
            crossed_junctions=[i for i, _, _ in crossed],
            junction_crossings=[dict(junction_index=i, boundaries=[kind for kind, _ in boundaries],
                                     termination_only=all(b >= coding_end for _, b in boundaries))
                                for i, _, boundaries in crossed],
            full_interval_support=dict(
                **rna_support(segments),
                read_lineage=lineage(segments) if lineage is not None else None,
                missing_quality_reads=len({observations[w["observation"]].identity for w in witnesses
                                           if observations[w["observation"]].missing_qualities}),
                scope="built_direct_junction_observations", witnesses=witnesses)))
    return dict(candidates=candidates, candidate_limit_reached=limited,
                scope="ATG_in_retained_path_crossing_event_related_junction",
                peptide_novelty_assessed=False)
