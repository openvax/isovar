"""Exploratory ATG ORFs and their observed RNA evidence, not initiation calls."""

from .genetic_code import standard_genetic_code
from .orf_inclusion import annotate_orf_inclusion
from .orf_start import annotate_orf_start
from .sv_rna_relations import EVENT_CROSSING_RELATIONS


def _witnesses(sequence, positions, start, end, junctions, observations):
    """Exact full-interval witnesses with compatible placements at this join.

    Sequence agreement includes unplaced bases. Two partial reads never
    become a full witness. Different assignments of junction homology are
    allowed only at unplaced bases, with shared placed anchors on each
    available side. A shorter ORF may need additional sequence to link it to
    the other flank; that interval must also match in this same observation.
    """
    dna = sequence[start:end]
    witnesses = {}
    for junction in junctions:
        left, right = junction["query_interval"]
        for key in junction["direct_observations"]:
            observation = observations[key]
            offset = observation.sequence.find(dna)
            while offset >= 0:
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
                    pairs = zip(positions[a:b], observation.positions[a + shift:b + shift])
                    consistent = (sequence[a:b] == observation.sequence[a + shift:b + shift]
                                  and all(p is None or q is None or p == q for p, q in pairs))
                else:
                    consistent = False
                if consistent:
                    link_interval = [a + shift, b + shift]
                    a = (observation.query_interval[1] - offset - len(dna) if observation.reverse
                         else observation.query_interval[0] + offset)
                    witness = witnesses.setdefault((key, offset), dict(
                        observation=key, observation_interval=[offset, offset + len(dna)],
                        original_query_interval=[a, a + len(dna)], junction_links=[]))
                    link = dict(junction_query_interval=[left, right], observation_interval=link_interval)
                    if link not in witness["junction_links"]:
                        witness["junction_links"].append(link)
                offset = observation.sequence.find(dna, offset + 1)
    return [witnesses[key] for key in sorted(witnesses)]


def _boundaries(junction):
    """Named RNA boundaries; unplaced bases need not be a DNA insertion."""
    left, right = junction["query_interval"]
    if right == left + 1 and junction["left"] is not None and junction["right"] is not None:
        return [("flank_to_flank", right)]
    return ([("donor_to_unplaced", left + 1)] if junction["left"] is not None else []) + (
        [("unplaced_to_acceptor", right)] if junction["right"] is not None else [])


def exploratory_orfs(sequence, positions, junctions, observations, references, cell_umi_support,
                     min_amino_acids, max_candidates, lineage=None, competing_splices=(), inclusion_thresholds=None):
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
        witnesses = _witnesses(sequence, positions, start, end, [j for _, j, _ in crossed], observations)
        segments = {observations[w["observation"]].identity for w in witnesses}
        labels = cell_umi_support(segments)
        start_evidence = annotate_orf_start(sequence, positions, start, references)
        start_evidence = annotate_orf_inclusion(
            start_evidence, sequence, positions, end, references, observations, witnesses,
            competing_splices=competing_splices, **(inclusion_thresholds or {}))
        for assessment in start_evidence["assessments"]:
            inclusion = assessment.get("splice_inclusion")
            if inclusion is not None:
                qualified = {observations[key].identity for key in inclusion["qualified_observations"]}
                inclusion["cell_umi_support"] = cell_umi_support(qualified)
                inclusion["read_lineage"] = lineage(qualified) if lineage is not None else None
        comparisons = start_evidence["reference_comparisons"]
        minus_three = sequence[start - 3] if start >= 3 else None
        plus_four = sequence[start + 3] if start + 3 < len(sequence) else None
        candidates.append(dict(
            query_interval=[start, end], amino_acids=aa, ends_with_stop_codon=has_stop,
            initiation_observed=False, translation_observed=False,
            start_codon_basis="observed_ATG", annotated_start=any(c["start_kind"] == "annotated_start" for c in comparisons),
            start_context=dict(sequence=sequence[max(0, start - 6):min(len(sequence), start + 9)],
                               query_interval=[max(0, start - 6), min(len(sequence), start + 9)],
                               minus_three=minus_three, plus_four=plus_four),
            reference_comparisons=comparisons, start_evidence=start_evidence,
            crossed_junctions=[i for i, _, _ in crossed],
            junction_crossings=[dict(junction_index=i, boundaries=[kind for kind, _ in boundaries],
                                     termination_only=all(b >= coding_end for _, b in boundaries))
                                for i, _, boundaries in crossed],
            full_interval_support=dict(
                segments=len(segments), fragments=len({s[:2] for s in segments}),
                molecule_labels=labels["complete_label_count"], cell_umi_support=labels,
                read_lineage=lineage(segments) if lineage is not None else None,
                missing_quality_segments=len({observations[w["observation"]].identity for w in witnesses
                                              if observations[w["observation"]].missing_qualities}),
                scope="built_direct_junction_observations", witnesses=witnesses)))
    return dict(candidates=candidates, candidate_limit_reached=limited,
                scope="ATG_in_retained_path_crossing_event_related_junction",
                peptide_novelty_assessed=False)
