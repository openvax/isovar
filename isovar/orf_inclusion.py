"""Observed splice linkage for intronic ATGs; no motif-based initiation calls."""

from copy import deepcopy

from .default_parameters import (
    SV_INCLUSION_MAPQ_255_IS_UNIQUE,
    SV_INCLUSION_MIN_BASE_QUALITY,
    SV_INCLUSION_MIN_MAPPING_QUALITY,
    SV_INCLUSION_MIN_SPLICE_ANCHOR_BASES,
)
from .orf_start import summarize_orf_start_evidence


def _runs(positions):
    runs = []
    for q, p in enumerate(positions):
        if p is None:
            continue
        previous = positions[q - 1] if q else None
        if (runs and runs[-1][1] == q and previous is not None
                and p == (previous[0], previous[1] + (1 if p[2] == "+" else -1), previous[2])):
            runs[-1][1] += 1
        else:
            runs.append([q, q + 1])
    return runs


def _edges(positions, reference):
    exons = list(reference.exons)
    if reference.strand == "-":
        exons.reverse()
    starts = [a if reference.strand == "+" else b - 1 for a, b in exons]
    ends = [b - 1 if reference.strand == "+" else a for a, b in exons]
    edges = []
    runs = _runs(positions)
    for (lo, a), (b, hi) in zip(runs, runs[1:]):
        left, right = positions[a - 1], positions[b]
        same = left[0] == right[0] == reference.contig and left[2] == right[2] == reference.strand
        forward = same and (right[1] - left[1]) * (1 if reference.strand == "+" else -1) > 0
        donor = ends.index(left[1]) if same and left[1] in ends else None
        acceptor = starts.index(right[1]) if same and right[1] in starts else None
        mechanism = "unresolved"
        if a == b and forward:
            if donor is not None and acceptor is not None:
                mechanism = "annotated_splice" if acceptor == donor + 1 else (
                    "exon_skipping" if acceptor > donor + 1 else "unresolved")
            elif donor is not None:
                mechanism = "cryptic_acceptor"
            elif acceptor is not None:
                mechanism = "cryptic_donor"
        elif not forward:
            mechanism = "rearranged_path"
        edges.append(dict(query_interval=[a - 1, b], left=list(left), right=list(right),
                          flanking_runs=[[lo, a], [b, hi]], mechanism=mechanism))
    return runs, edges


def _qualified_link(sequence, positions, start, end, edge, witness, observation, min_anchor, min_baseq, min_mapq,
                    mapq_255_is_unique):
    a, b = edge["query_interval"]
    left_run, right_run = edge["flanking_runs"]
    if a - left_run[0] + 1 < min_anchor or right_run[1] - b < min_anchor:
        return None
    shift = witness["observation_interval"][0] - start
    # Both the ORF's event link and the transcript-inclusion splice must be
    # carried by this same observation, with exact intervening sequence.
    links = witness.get("junction_links", [])
    required = edge.get("required_interval", [start, end])
    lo = min(start, required[0], a + 1 - min_anchor,
             *(x["observation_interval"][0] - shift for x in links))
    hi = max(end, required[1], b + min_anchor,
             *(x["observation_interval"][1] - shift for x in links))
    x, y = lo + shift, hi + shift
    row = dict(observation=witness["observation"], query_interval=[lo, hi], observation_interval=[x, y],
               minimum_mapping_quality=observation.minimum_mapping_quality,
               mapping_quality_255=observation.mapping_quality_255,
               minimum_base_quality=None, status="unlinked")
    if (x < 0 or y > len(observation.sequence) or sequence[lo:hi] != observation.sequence[x:y]
            or any(p is not None and p != q for p, q in zip(positions[lo:hi], observation.positions[x:y]))):
        return row
    observed_kind = next((kind for i, j, kind in observation.breaks if (i, j) == (a + shift, b + shift)), None)
    row["alignment_gap_kind"] = observed_kind
    if observed_kind != "N":
        row["status"] = "no_observed_splice"
        return row
    qualities = observation.qualities[x:y]
    if len(qualities) == y - x and all(q is not None for q in qualities):
        row["minimum_base_quality"] = min(qualities)
    mapq = observation.minimum_mapping_quality
    row["status"] = (
        "mapping_quality_unavailable" if observation.mapping_quality_255 and not mapq_255_is_unique else
        "low_mapping_quality" if mapq is not None and mapq < min_mapq else
        "base_quality_unavailable" if row["minimum_base_quality"] is None else
        "low_base_quality" if row["minimum_base_quality"] < min_baseq else "qualified")
    return row


def annotate_orf_inclusion(annotation, sequence, positions, end, references, observations, witnesses,
                           *, min_splice_anchor_bases=SV_INCLUSION_MIN_SPLICE_ANCHOR_BASES,
                           min_base_quality=SV_INCLUSION_MIN_BASE_QUALITY,
                           min_mapping_quality=SV_INCLUSION_MIN_MAPPING_QUALITY,
                           mapq_255_is_unique=SV_INCLUSION_MAPQ_255_IS_UNIQUE,
                           competing_splices=()):
    """Link an intronic ATG, observed splice and SV on the same RNA observation.

    Parameters
    ----------
    annotation : dict
        Result of ``annotate_orf_start``; copied without modifying the input.
    sequence, positions : str, sequence
        Retained path and its per-base placements in candidate orientation.
    end : int
        ORF end, including a stop codon if observed.
    references : iterable of FusionReference
        The same versioned transcript models used for start annotation.
    observations : mapping
        Keys to ``sv_rna.RnaObservation`` objects, including original-record
        qualities, mapping quality and CIGAR-derived gap kinds.
    witnesses : iterable of dict
        Full-ORF, event-linked witnesses from RNA reconstruction. Partial
        observations cannot substitute for same-observation linkage.
    min_splice_anchor_bases, min_base_quality, min_mapping_quality : int
        Explicit conservative thresholds for transcript-inclusion promotion:
        exact bases on each side of the splice, and the minimum base and
        mapping quality over the linked interval.
    mapq_255_is_unique : bool
        Whether MAPQ 255 means a unique alignment, as STAR writes it (the
        default), rather than the SAM specification's "unavailable", which
        fails the gate.
    competing_splices : iterable of dict
        Observed CIGAR N joins with ``left``, ``right``, and ``fragments`` in
        this input source, after MAPQ filtering. Missing competitors mean
        retention is unresolved, not absent. Counts are never summed here.

    Returns
    -------
    dict
        Transcript-specific path graph and qualified/unqualified linkage.
        Tier 3 requires a CIGAR N join from the ATG's contiguous intronic
        block to an annotated exon boundary, exact linked sequence, and
        known sufficient mapping/base quality. Two such joins describe
        exonization. Intronic coverage, deletion gaps, supplementary joins,
        motifs and unrelated splices never promote an ATG. Retention also
        requires both exon/intron boundaries, a separate annotated splice
        on the same observation and a competing spliced path. Rearranged
        intronic sequence requires a linked annotated partner exon and
        an annotated splice on that same observation. This does not infer RNA strand,
        initiation, translation, or independence of templates.
    """
    for value in (min_splice_anchor_bases, min_base_quality, min_mapping_quality):
        if type(value) is not int or value < 1:
            raise ValueError("Inclusion thresholds must be positive integers")
    if not isinstance(mapq_255_is_unique, bool):
        raise ValueError("mapq_255_is_unique must be True or False")
    start = annotation["query_offset"]
    if len(sequence) != len(positions) or not 0 <= start < end <= len(sequence):
        raise ValueError("Invalid ORF/path intervals")
    result = deepcopy(annotation)
    references = {(r.transcript_id, r.annotation): r for r in references}
    witnesses = list(witnesses)
    competing_splices = list(competing_splices)
    for assessment in result["assessments"]:
        reference = references.get((assessment["transcript_id"], assessment["annotation"]))
        # Only an intronic start needs splice evidence for inclusion; other
        # regions keep splice_inference="not_assessed" and no inclusion record.
        if reference is None or assessment["region"] != "intronic":
            continue
        runs, edges = _edges(positions, reference)
        evidence = dict(policy="isovar.orf_inclusion.v1", status="unresolved", mechanisms=[],
                        graph=dict(reference_exons=[list(e) for e in reference.exons],
                                   observed_runs=runs, observed_edges=edges),
                        thresholds=dict(min_splice_anchor_bases=min_splice_anchor_bases,
                                        min_base_quality=min_base_quality,
                                        min_mapping_quality=min_mapping_quality,
                                        mapq_255_is_unique=mapq_255_is_unique),
                        witnesses=[], qualified_observations=[], qualified_fragments=0,
                        competing_splices=[],
                        matched_normal="not_assessed", splice_prediction="not_assessed",
                        mature_transcript_proven=False)
        assessment["splice_inclusion"] = evidence
        assessment["splice_inference"] = "observed_path_assessed"
        lo, hi = assessment["intron_interval"]
        start_run = next((r for r in runs if r[0] <= start and start + 3 <= r[1]), None)
        relevant = [e for e in edges if
                    (e["mechanism"] == "cryptic_acceptor" and e["flanking_runs"][1] == start_run
                     and lo <= e["right"][1] < hi) or
                    (e["mechanism"] == "cryptic_donor" and e["flanking_runs"][0] == start_run
                     and lo <= e["left"][1] < hi)]
        # Coverage becomes a retention hypothesis only when a single aligned
        # run crosses both boundaries. Processed context and a competitor
        # are additional requirements, not substitutes for those boundaries.
        if start_run is not None:
            before, after = ((lo - 1, hi) if reference.strand == "+" else (hi, lo - 1))
            offset = {p: q for q, p in enumerate(positions[start_run[0]:start_run[1]], start_run[0])}
            left = (reference.contig, before, reference.strand)
            right = (reference.contig, after, reference.strand)
            a, b = offset.get(left), offset.get(right)
            if (a is not None and b is not None and a - start_run[0] + 1 >= min_splice_anchor_bases
                    and start_run[1] - b >= min_splice_anchor_bases):
                competitors = [c for c in competing_splices if c["left"] == list(left)
                               and c["right"] == list(right) and c["fragments"] > 0]
                evidence["competing_splices"] = deepcopy(competitors)
                if competitors:
                    relevant.extend(dict(e, mechanism="intron_retention",
                                         required_interval=[a + 1 - min_splice_anchor_bases, b + min_splice_anchor_bases])
                                    for e in edges if e["mechanism"] == "annotated_splice")
                else:
                    evidence["status"] = "retention_without_competing_splice"
            # A rearrangement alone is not mature-transcript evidence. Demand
            # a linked partner exon with its own observed annotated splice.
            for edge in edges:
                if edge["mechanism"] != "rearranged_path" or start_run not in edge["flanking_runs"]:
                    continue
                own = edge["left"] if edge["flanking_runs"][0] == start_run else edge["right"]
                if not lo <= own[1] < hi:
                    continue
                other = edge["right"] if edge["flanking_runs"][0] == start_run else edge["left"]
                for partner in references.values():
                    if (other[0] != partner.contig or other[2] != partner.strand
                            or partner.spliced_offset(other[1]) is None):
                        continue
                    _, context = _edges(positions, partner)
                    required = [edge["query_interval"][0] + 1 - min_splice_anchor_bases,
                                edge["query_interval"][1] + min_splice_anchor_bases]
                    if required[0] < edge["flanking_runs"][0][0] or required[1] > edge["flanking_runs"][1][1]:
                        continue
                    relevant.extend(dict(e, mechanism="fusion_rearranged_exon", required_interval=required,
                                         context_transcript=partner.transcript_id)
                                    for e in context if e["mechanism"] == "annotated_splice")
        qualified = set()
        for edge in relevant:
            for witness in witnesses:
                obs = observations[witness["observation"]]
                row = _qualified_link(sequence, positions, start, end, edge, witness, obs,
                                      min_splice_anchor_bases, min_base_quality, min_mapping_quality,
                                      mapq_255_is_unique)
                if row is None:
                    continue
                row["mechanism"] = edge["mechanism"]
                if "context_transcript" in edge:
                    row["context_transcript"] = edge["context_transcript"]
                row["splice_query_interval"] = edge["query_interval"]
                if row not in evidence["witnesses"]:
                    evidence["witnesses"].append(row)
                if row["status"] == "qualified":
                    qualified.add(witness["observation"])
        evidence["qualified_observations"] = sorted(qualified)
        evidence["qualified_fragments"] = len({observations[key].fragment for key in qualified})
        mechanisms = {w["mechanism"] for w in evidence["witnesses"] if w["status"] == "qualified"}
        evidence["mechanisms"] = sorted(mechanisms)
        if mechanisms == {"cryptic_acceptor", "cryptic_donor"} and any(
                {w["mechanism"] for w in evidence["witnesses"]
                 if w["observation"] == key and w["status"] == "qualified"} == mechanisms for key in qualified):
            evidence["mechanisms"].append("exonization")
        if qualified:
            evidence["status"] = "same_observation_splice_linked"
            assessment.update(tier="intronic_splice_supported", priority=3,
                              transcript_inclusion="same_observation_splice_linked")
    result.update(summarize_orf_start_evidence([result]))
    return result
