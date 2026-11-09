"""Shared exact RNA path discovery and complete-input original-fragment support.

Regions are contig -> (zero-based start, exclusive end, orientation) intervals.
Reference transcript models and coding interpretation belong to callers.
"""

from collections import Counter, defaultdict
from dataclasses import replace

from .cell_umi import CellUmiEvidence
from .chimeric_alignment import (
    alignment_placement,
    matches_reported_alignment,
    sa_entries,
)
from .read_identity import sam_record_id, segment_identity
from .read_lineage import ReadLineage
from .read_metadata import record_evidence
from .rna_observation import (
    ObservationGeometryCache,
    ObservationIndex,
    RnaObservation,
    build_segment_observations,
)
from .rna_path_assembly import extend_path_both_directions
from .rna_path_geometry import (
    PathPositions,
    PlacementIndex,
    exact_overlap_agrees,
    overlap_offsets,
    position_runs,
    reverse_path,
)
from .rna_path_linkage import OriginalPathMatch, template_path_linkage


def matching_observation_intervals(sequence, positions, observation, *, contained=False):
    """Find exact, coordinate-consistent overlaps with an original observation.

    Parameters
    ----------
    sequence : str
        Observed RNA bases in the orientation of positions; no reference filling.
    positions : sequence of Position or None
        One placement per query base: (contig, zero-based base, strand), or
        None for an observed base whose placement is unresolved.
    observation : RnaObservation
        One original segment and placement choice, retaining its identity.
    contained : bool, optional
        Require the entire observation inside the candidate path (default False).

    Returns
    -------
    iterator of tuple[int, int]
        Half-open intervals on the candidate path. Shared placed anchors route
        the match; exact bases and every jointly placed coordinate must agree.
    """
    for _, offset in sorted(overlap_offsets(position_runs(positions), 0, observation.key, position_runs(observation.positions))):
        if contained and (offset < 0 or offset + len(observation.sequence) > len(sequence)):
            continue
        lo, hi = max(0, offset), min(len(sequence), offset + len(observation.sequence))
        if lo < hi and exact_overlap_agrees(sequence, positions, observation, offset, lo, hi):
            yield lo, hi


def path_contains(first, second):
    """Check exact sequence and placement containment of two RNA paths.

    Parameters
    ----------
    first, second : tuple[str, sequence]
        (Observed sequence, per-base placements) paths in the same orientation.

    Returns
    -------
    bool
        Whether second is wholly contained in first with identical placements.
    """
    sequence, positions = first
    other, mapping = second
    start = sequence.find(other)
    while start >= 0:
        if positions[start:start + len(other)] == mapping:
            return True
        start = sequence.find(other, start + 1)
    return False


def _relevant(read, regions):
    return not read.is_unmapped and read.reference_end is not None and any(
        read.reference_start < end and read.reference_end > start
        for start, end, _ in regions.get(read.reference_name, ()))


def _target_orientation(observation, regions):
    return any(p is not None and any(start <= p[1] < end and p[2] == strand
                                    for start, end, strand in regions.get(p[0], ()))
               for p in observation.positions)


def original_template_observations(template, collector, reasons, regions, cache):
    """Build target-oriented alternatives for the original segments of a template.

    Parameters
    ----------
    template : sequence of pysam.AlignedSegment
        Eligible original records from the library index, including competing
        alignments and mates. Records remain grouped by original segment.
    collector : ReadCollector
        Read/end processing policy used to prepare the library index.
    reasons : mutable counter
        Receives observation rejection and limitation counts.
    regions : mapping
        Contig to zero-based, half-open (start, end, orientation) intervals.
    cache : ObservationGeometryCache or None
        Optional bounded geometry cache; placement interning is per template.

    Yields
    ------
    tuple
        Original segment identity and its target-oriented observations. Segments
        with an overlapping alignment are yielded even if all their observation
        alternatives are rejected, in which case the observation list is empty.
    """
    grouped = defaultdict(list)
    for read in template:
        grouped[segment_identity(read)].append(read)
    for identity, records in sorted(grouped.items()):
        if any(_relevant(read, regions) for read in records):
            # The placement interner is per template: full-input support does
            # not retain every background genomic base in an unbounded cache.
            yield identity, [o for o in build_segment_observations(identity, records, collector, reasons, {}, cache)
                             if _target_orientation(o, regions)]


def _seed_choices(observations, candidates=None):
    """Share geometric seeds only when an unambiguous original can seed them.

    All observations remain in the overlap index. An unambiguous parent has no
    other placement that could extend its own seed, so its identical geometric
    copies are interchangeable here. Ambiguous parents retain their separate
    seed constraints. Reciprocal orientation views are one placement choice.
    """
    placements, geometries = defaultdict(set), defaultdict(list)
    for observation in observations:
        placements[observation.identity].add(observation.key[:-1])
    for observation in observations if candidates is None else candidates:
        geometries[observation.comparison_key].append(observation)
    seeds = []
    for members in geometries.values():
        unambiguous = [o for o in members if len(placements[o.identity]) == 1]
        seeds.extend([min(unambiguous, key=lambda o: o.key)] if unambiguous else members)
    return sorted(seeds, key=lambda o: (-len(o.sequence), o.sequence, repr(o.positions), o.key))


def assemble_rna_paths(chosen, max_paths, min_overlap, *, seed_observations: tuple[RnaObservation, ...] | None = None):
    """Assemble maximal exact paths from original observation alternatives.

    Parameters
    ----------
    chosen : mapping
        Original observation keys to RnaObservation objects.
    max_paths : int
        Bound on retained discovered paths; reaching it is reported.
    min_overlap : int
        Minimum exact overlap for extension.
    seed_observations : tuple of RnaObservation, optional
        Explicit event-adapter seeds. Their original placement constraints
        remain attached through both extension directions.

    Returns
    -------
    tuple[list, set[str]]
        Maximal (sequence, placement_tuple) paths and limit/ambiguity notes.

    Notes
    -----
    This is discovery, not quantification. Identical geometry can share work,
    but alternative placements of one parent cannot be independent witnesses.
    """
    observations = {o.key: o for o in chosen.values()}
    if seed_observations is not None:
        observations.update((o.key, o) for o in seed_observations)
    originals = tuple(observations.values())
    # extend_path_both_directions mirrors the path to extend its 5' end. The index needs the
    # same mirror views, retaining the original segment/fragment identity.
    for o in originals:
        sequence, positions = reverse_path(o.sequence, o.positions)
        key = o.key[:-1] + ("-" if o.key.endswith("+") else "+")
        n = len(o.sequence)
        observations.setdefault(key, replace(
            o, key=key, sequence=sequence, positions=positions,
            breaks=tuple(sorted((n - 1 - j, n - 1 - i, kind) for i, j, kind in o.breaks)),
            reverse_complement=not o.reverse_complement, qualities=o.qualities[::-1]))
    index = ObservationIndex(list(observations.values()))
    paths, notes = {}, set()
    for seed in _seed_choices(originals, seed_observations):
        key = seed.sequence, tuple(seed.positions)
        if len(paths) >= max_paths:
            notes.add("path_limit")
            break
        for sequence, positions, _, _ in extend_path_both_directions(
                *key, index, min_overlap, (0, 0, 0), max_paths - len(paths), [], notes,
                seed_observation=seed):
            paths[sequence, tuple(positions)] = None
            if len(paths) >= max_paths:
                notes.add("path_limit")
                break
    maximal = [p for p in paths if not any(len(q[0]) > len(p[0]) and path_contains(q, p) for q in paths)]
    return maximal, notes


def discover_rna_paths(groups, regions, collector, max_discovery_reads, max_paths, min_overlap, *,
              seed_observations: tuple[RnaObservation, ...] | None = None):
    """Discover bounded observed RNA paths from an original-library index.

    Parameters
    ----------
    groups : RnaLibraryIndex
        Open original-library index with matching source and eligibility policy.
    regions : mapping
        Contig to (start, end, strand) search intervals; bounds are zero-based
        and half-open. Strand chooses path orientation, not biological RNA strand.
    collector : ReadCollector
        Explicit read/end policy, applied to original alignments.
    max_discovery_reads : int
        Maximum eligible original segments selected for assembly, not a support cap.
    max_paths : int
        Bound on distinct discovered paths.
    min_overlap : int
        Minimum exact overlap used by assemble_rna_paths.
    seed_observations : tuple of RnaObservation, optional
        Adapter-selected original seeds, retaining their placement constraints.

    Returns
    -------
    tuple[list, dict]
        Observed (sequence, placement_tuple) paths and discovery accounting:
        considered/eligible/selected segments, budgets, rejections and truncation.

    Notes
    -----
    The caller owns the open library index. A truncated search does not prove
    absence of another path, nor identify a dominant isoform.
    """
    chosen, reasons, cache = {}, Counter(), ObservationGeometryCache()
    considered, eligible, selected = 0, 0, 0
    for template in groups.templates([(contig, start, end) for contig, spans in regions.items()
                                      for start, end, _ in spans]):
        for _, observations in original_template_observations(template, collector, reasons, regions, cache):
            considered += 1
            if not observations:
                continue
            eligible += 1
            if selected >= max_discovery_reads:
                continue
            selected += 1
            for o in observations:
                # Geometry can share matching work, but it cannot substitute
                # for the original parent when enforcing placement choices.
                chosen.setdefault(o.key, o)
    maximal, notes = assemble_rna_paths(chosen, max_paths, min_overlap, seed_observations=seed_observations)
    return maximal, dict(template_order="read_group_query_name", target_segments=considered,
                         eligible_observed_segments=eligible, selected_segments=selected,
                         distinct_observations=len({o.comparison_key for o in chosen.values()}),
                         max_discovery_reads=max_discovery_reads,
                         max_paths=max_paths, min_overlap=min_overlap,
                         # Truncation means a budget was reached; other notes
                         # (e.g. retained alternative placements) are limitations only.
                         truncated=eligible > selected or "path_limit" in notes, limitations=sorted(notes),
                         observation_rejections=dict(reasons))


def support_rna_paths(groups, candidates, regions, collector, header, sample_id, source, *, label="isoform"):
    """Measure fixed candidate paths against the complete declared original input.

    Parameters
    ----------
    groups : RnaLibraryIndex
        Open original-library index with matching source and eligibility policy.
    candidates : list of dict
        Mutable native candidate records: id, sequence, immutable placements,
        and coding entries. Complete coding entries need status=complete_CDS
        and half-open query_start/query_end. Optional features enable linkage.
    regions : mapping
        Contig to (start, end, strand) search intervals; bounds are zero-based
        and half-open. Strand chooses path orientation, not biological RNA strand.
    collector : ReadCollector
        Explicit read/end policy, applied to original alignments.
    header : dict
        Original SAM header, preserving read groups and processing lineage.
    sample_id, source : str
        Original evidence scope; counts from different sources are not pooled.
    label : str, optional
        Existing native schema label (isoform, sequence or path). Changes
        identifier-key spelling only, never attribution or counting rules.

    Returns
    -------
    dict
        Complete-input target/unassigned fragment totals, equivalence classes,
        original witnesses and observation rejections. candidates are also updated
        with support fields; no new sidecar or alternate schema is created.

    Notes
    -----
    supporting_segments/fragments and rna_support count partial-overlap
    contributors. coverage is per-base original-segment depth. compatible_fragments
    requires compatible observations for a template; unique_fragments additionally
    requires one candidate and no unresolved outside placement. Whole-path/CDS
    segment witnesses are separate from tiled coverage. Shared alternative counts
    are not additive. Discovery truncation does not cap this support pass.
    """
    compatible_key, overlapping_key = "compatible_" + label + "s", "overlapping_" + label + "s"
    ids_key = label + "_ids"
    coverage = {c["id"]: [0] * len(c["sequence"]) for c in candidates}
    by_id = {c['id']: c for c in candidates}
    # Reuse each immutable candidate's collinear blocks across the entire
    # original-input recount. Raw output placements remain unchanged.
    placements = {c["id"]: PathPositions(c["placements"]) for c in candidates}
    candidate_index = PlacementIndex({i: c['placements'] for i, c in enumerate(candidates)})
    supporters, whole, assignments = defaultdict(set), defaultdict(set), Counter()
    whole_cds = defaultdict(set)
    linkage = {c["id"]: [] for c in candidates if "features" in c}
    witnesses, reasons, cache = [], Counter(), ObservationGeometryCache()
    lineage = ReadLineage(groups, header)
    labels = CellUmiEvidence(groups, header, sample_id, source, template_loader=groups.template)
    target_fragments, unassigned = 0, 0
    for template in groups.templates([(contig, start, end) for contig, spans in regions.items()
                                      for start, end, _ in spans]):
        observations = list(original_template_observations(template, collector, reasons, regions, cache))
        if not observations:
            continue
        target_fragments += 1
        fragment_matches, segment_rows = None, []
        template_matches = defaultdict(list)
        outside = any(not _relevant(read, regions) for read in template)
        uncertainty = placement_uncertainty(template)
        outside |= bool(uncertainty)
        for identity, observed in observations:
            compatible, supporting = set(), set()
            # An exact overlap needs a shared placed base. Use the same run
            # index as assembly, then retain the original exact checks below.
            relevant = sorted({i for o in observed for i, _ in
                               candidate_index.overlapping(o.positions, 0, len(o.positions))})
            for i in relevant:
                candidate = candidates[i]
                cid = candidate["id"]
                covered, full = set(), False
                for o in observed:
                    for lo, hi in matching_observation_intervals(candidate["sequence"], placements[cid], o):
                        if cid in linkage:
                            template_matches[cid].append(OriginalPathMatch(identity, o.key, lo, hi))
                        covered.update(range(lo, hi))
                        full |= lo == 0 and hi == len(candidate["sequence"])
                        for i, coding in enumerate(candidate["coding"]):
                            if coding["status"] == "complete_CDS" and lo <= coding["query_start"] and hi >= coding["query_end"]:
                                whole_cds[cid, i].add(identity)
                    if any(matching_observation_intervals(candidate["sequence"], placements[cid], o, contained=True)):
                        compatible.add(cid)
                if covered:
                    supporting.add(cid)
                    supporters[cid].add(identity)
                    for offset in covered:
                        coverage[cid][offset] += 1
                if full:
                    whole[cid].add(identity)
            fragment_matches = compatible if fragment_matches is None else fragment_matches & compatible
            segment_rows.append(dict(identity=list(identity), **{compatible_key: sorted(compatible), overlapping_key: sorted(supporting)},
                                     lineage=lineage.segment(identity), cell_umi=labels.segment(identity)))
        ids = tuple(sorted(fragment_matches or ()))
        for cid, matches in template_matches.items():
            if matches:
                candidate = by_id[cid]
                linkage[cid].append(dict(fragment=list(observations[0][0][:2]),
                    **template_path_linkage(matches, len(candidate["sequence"]), candidate["features"])))
        assignments[ids, outside] += 1
        unassigned += not ids
        witnesses.append(dict(fragment=list(observations[0][0][:2]), **{compatible_key: list(ids)},
                              competing_placement_outside_targets=outside, placement_uncertainty=uncertainty, segments=segment_rows,
                              records=[dict(id=sam_record_id(read), contig=read.reference_name,
                                            start=read.reference_start, cigar=read.cigarstring,
                                            supplementary_alignments=read.get_tag("SA") if read.has_tag("SA") else None,
                                            mate_contig=read.next_reference_name if read.is_paired and read.next_reference_id >= 0 else None,
                                            mate_start=read.next_reference_start if read.is_paired else None,
                                            **record_evidence(read)) for read in template]))
        # Witnesses hold the consulted metadata. Caches need not grow with the
        # total file; ancestry can still consult the disk-backed source index.
        lineage.cache.clear()
        labels.cache.clear()
    for candidate in candidates:
        cid = candidate["id"]
        candidate.update(coverage=coverage[cid], supporting_segments=len(supporters[cid]),
                         supporting_fragments=len({key[:2] for key in supporters[cid]}),
                         whole_path_segments=[list(key) for key in sorted(whole[cid])],
                         unique_fragments=assignments[((cid,), False)],
                         compatible_fragments=sum(n for (ids, _), n in assignments.items() if cid in ids))
        candidate["sequence_status"] = (
            "whole_path_observed" if whole[cid] else "fully_covered_hypothesis" if all(coverage[cid])
            else "partially_covered_reference" if any(coverage[cid]) else "not_observed")
        for i, coding in enumerate(candidate["coding"]):
            coding["whole_cds_segments"] = [list(key) for key in sorted(whole_cds[cid, i])]
            coding["all_bases_observed"] = all(coverage[cid][coding["query_start"]:coding["query_end"]])
        candidate["rna_support"] = labels.support(supporters[cid])
        candidate["signal_support"] = lineage.support(supporters[cid])
        if cid in linkage:
            rows = linkage[cid]
            candidate["linkage"] = dict(
                single_read_whole_path_fragments=sum(bool(row["whole_path_single_read_segments"]) for row in rows),
                paired_template_whole_path_fragments=sum(row["whole_path_paired_template"] for row in rows),
                paired_template_only_whole_path_fragments=sum(
                    row["whole_path_paired_template"] and not row["whole_path_single_read_segments"] for row in rows),
                witnesses=rows)
        labels.cache.clear()
        lineage.cache.clear()
    classes = [dict(**{ids_key: list(ids)}, outside_targets=outside, fragments=n)
               for (ids, outside), n in sorted(assignments.items())]
    return dict(scope="complete_declared_input", target_fragments=target_fragments,
                unassigned_target_fragments=unassigned, equivalence_classes=classes,
                observation_rejections=dict(reasons), witnesses=witnesses)


def placement_uncertainty(template):
    """Describe unavailable or uncertain placements preventing unique attribution.

    Parameters
    ----------
    template : sequence of pysam.AlignedSegment
        Actual original records for one template, including supplied alternative
        alignments and mates. Tags alone never supply missing observations.

    Returns
    -------
    list of str
        Sorted reason codes for unavailable MAPQ, missing mapped mates or SA/NH
        placements, zero MAPQ, or malformed SA tags. An empty list does not by
        itself establish unique support; candidate compatibility is separate.
    """
    reasons, segments = set(), defaultdict(list)
    for read in template:
        segments[segment_identity(read)].append(read)
    for identity, records in segments.items():
        for read in records:
            if read.mapping_quality == 0:
                reasons.add("mapping_quality_zero")
            if read.mapping_quality == 255 and not read.has_tag("NH"):
                reasons.add("mapping_quality_unavailable")
            if read.has_tag("NH") and read.get_tag("NH") > sum(not r.is_supplementary for r in records):
                reasons.add("unreported_NH_placements")
            if read.is_paired and not read.mate_is_unmapped:
                mate = 128 if read.is_read1 else 64
                if (*identity[:2], mate) not in segments:
                    reasons.add("missing_mapped_mate")
            if read.has_tag("SA"):
                try:
                    for contig, start, strand, cigar, mapq, nm in sa_entries(read.get_tag("SA")):
                        if strand not in {"+", "-"} or not 0 <= mapq <= 255 or nm < 0:
                            raise ValueError("Invalid SA entry")
                        reported = alignment_placement(contig, start, strand == "-", cigar)
                        if not any(r is not read and not r.is_unmapped and matches_reported_alignment(
                                alignment_placement(r.reference_name, r.reference_start, r.is_reverse, r.cigarstring), reported)
                                   for r in records):
                            reasons.add("missing_SA_placement")
                except ValueError:
                    reasons.add("unparseable_SA_tag")
    return sorted(reasons)
