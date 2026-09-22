"""Exploratory ATG ORFs and their observed RNA evidence, not initiation calls."""

from .genetic_code import standard_genetic_code
from .read_metadata import record_evidence as record_evidence


def _start_references(sequence, positions, start, amino_acids, stop, models):
    comparisons = []
    for model in models:
        reference = model.reference
        offset = model.offsets.get(positions[start])
        if (offset is None or model.positions[offset:offset + 3] != list(positions[start:start + 3])
                or reference.sequence[offset:offset + 3] != "ATG"):
            continue
        if reference.cds_start is None:
            kind = "noncoding_transcript"
        elif offset == reference.cds_start:
            kind = "annotated_start"
        elif offset < reference.cds_start:
            kind = "five_prime_UTR"
        elif offset >= reference.cds_end:
            kind = "three_prime_UTR"
        else:
            kind = "internal_CDS"
        ref_aa, ref_stop = standard_genetic_code.translate(reference.sequence[offset:], first_codon_is_start=True)
        shared = 0
        for a, b in zip(amino_acids, ref_aa):
            if a != b:
                break
            shared += 1
        comparisons.append(dict(
            transcript_id=reference.transcript_id, annotation=reference.annotation,
            transcript_offset=offset, start_kind=kind, annotated_cds_start=reference.cds_start,
            amino_acids=ref_aa, ends_with_stop_codon=ref_stop,
            shared_prefix_amino_acids=shared,
            differs_from_reference_orf=(amino_acids, stop) != (ref_aa, ref_stop)))
    return comparisons


def _witnesses(sequence, positions, start, end, junctions, observations):
    """Exact full-interval witnesses with compatible placements at this join.

    Sequence agreement includes unplaced bases. Two partial reads never
    become a full witness. Different assignments of junction homology are
    allowed only at unplaced bases, with shared placed anchors on each
    available side. The observation already makes the cited junction.
    """
    dna = sequence[start:end]
    witnesses = {}
    for junction in junctions:
        left, right = junction["query_interval"]
        for key in junction["direct_observations"]:
            observation = observations[key]
            offset = observation.sequence.find(dna)
            while offset >= 0:
                pairs = list(zip(positions[start:end], observation.positions[offset:offset + len(dna)]))
                shared = [start + i for i, (p, q) in enumerate(pairs) if p is not None and p == q]
                consistent = all(p is None or q is None or p == q for p, q in pairs)
                sides = (any(q <= left for q in shared) if junction["left"] is not None else True,
                         any(q >= right for q in shared) if junction["right"] is not None else True)
                if consistent and shared and all(sides):
                    a = (observation.query_interval[1] - offset - len(dna) if observation.reverse
                         else observation.query_interval[0] + offset)
                    witnesses[key, offset] = dict(observation=key, observation_interval=[offset, offset + len(dna)],
                                                  original_query_interval=[a, a + len(dna)])
                offset = observation.sequence.find(dna, offset + 1)
    return [witnesses[key] for key in sorted(witnesses)]


def exploratory_orfs(sequence, positions, junctions, observations, models, cell_umi_support,
                     min_amino_acids, max_candidates, lineage=None):
    """Enumerate bounded ATG candidates crossing an event-related RNA join.

    These are potential translations of observed sequence, separate from
    annotated CDS frame transfer. Context is reported, never scored as proof
    of initiation. A missing stop means the retained sequence ends first.
    """
    related = [j for j in junctions if not j["annotated"] and j["relation"] in (
        "breakpoint_junction", "event_compatible_junction", "breakpoint_clip_partner_unplaced",
        "splice_ambiguous_event_junction")]
    candidates, limited = [], False
    # Next in-frame stop for each ATG, without translating every long suffix.
    next_stop = [None, None, None]
    starts = []
    for q in range(len(sequence) - 3, -1, -1):
        codon = sequence[q:q + 3]
        if codon in ("TAA", "TAG", "TGA"):
            next_stop[q % 3] = q
        elif codon == "ATG":
            stop = next_stop[q % 3]
            coding_end = stop if stop is not None else q + 3 * ((len(sequence) - q) // 3)
            crossed = [j for j in related if q <= j["query_interval"][0]
                       and j["query_interval"][1] < coding_end]
            if crossed and (coding_end - q) // 3 >= min_amino_acids:
                starts.append((q, coding_end, stop is not None, crossed))
    limited = len(starts) > max_candidates
    for start, coding_end, has_stop, crossed in reversed(starts[-max_candidates:]):
        end = coding_end + (3 if has_stop else 0)
        aa, _ = standard_genetic_code.translate(sequence[start:coding_end], first_codon_is_start=True)
        witnesses = _witnesses(sequence, positions, start, end, crossed, observations)
        segments = {observations[w["observation"]].identity for w in witnesses}
        labels = cell_umi_support(segments)
        comparisons = _start_references(sequence, positions, start, aa, has_stop, models)
        minus_three = sequence[start - 3] if start >= 3 else None
        plus_four = sequence[start + 3] if start + 3 < len(sequence) else None
        candidates.append(dict(
            query_interval=[start, end], amino_acids=aa, ends_with_stop_codon=has_stop,
            initiation_observed=False, translation_observed=False,
            start_codon_basis="observed_ATG", annotated_start=any(c["start_kind"] == "annotated_start" for c in comparisons),
            start_context=dict(sequence=sequence[max(0, start - 6):min(len(sequence), start + 9)],
                               query_interval=[max(0, start - 6), min(len(sequence), start + 9)],
                               minus_three=minus_three, plus_four=plus_four),
            reference_comparisons=comparisons,
            crossed_junctions=[junctions.index(j) for j in crossed],
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
