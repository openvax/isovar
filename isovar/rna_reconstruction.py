"""Shared exact RNA path discovery and complete-input original-fragment support.

Regions are contig -> (zero-based start, exclusive end, orientation) intervals.
Reference transcript models and coding interpretation belong to callers.
"""

from collections import Counter, defaultdict
from dataclasses import replace

from .cell_umi import CellUmiEvidence
from .chimeric_alignment import _alignment, _matches_reported, sa_entries
from .read_identity import segment_identity
from .read_lineage import ReadLineage
from .read_metadata import record_evidence
from .sv_rna import (
    _GeometryCache, _ObservationIndex, _agrees, _build_segment, _extend_both,
    _mirror, _offsets, _record_id, _runs,
)



def _matches(sequence, positions, observation, *, contained=False):
    """All exact, coordinate-consistent overlaps with an original observation."""
    for _, offset in sorted(_offsets(_runs(positions), 0, observation.key, _runs(observation.positions))):
        if contained and (offset < 0 or offset + len(observation.sequence) > len(sequence)):
            continue
        lo, hi = max(0, offset), min(len(sequence), offset + len(observation.sequence))
        if lo < hi and _agrees(sequence, positions, observation, offset, lo, hi):
            yield lo, hi


def _contains(first, second):
    """Whether an entire (sequence, placements) is contained exactly in another."""
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


def _observations(template, collector, reasons, regions, cache):
    grouped = defaultdict(list)
    for read in template:
        grouped[segment_identity(read)].append(read)
    for identity, records in sorted(grouped.items()):
        if any(_relevant(read, regions) for read in records):
            # The placement interner is per template: full-input support does
            # not retain every background genomic base in an unbounded cache.
            yield identity, [o for o in _build_segment(identity, records, collector, reasons, {}, cache)
                             if _target_orientation(o, regions)]


def _assemble(chosen, max_paths, min_overlap):
    """Maximal exact paths from distinct coordinate-aware observations."""
    observations = {o.key: o for o in chosen.values()}
    # _extend_both mirrors the path to extend its 5' end. The index needs the
    # same mirror views, retaining the original segment/fragment identity.
    for o in chosen.values():
        sequence, positions = _mirror(o.sequence, o.positions)
        key = o.key[:-1] + ("-" if o.key.endswith("+") else "+")
        n = len(o.sequence)
        observations.setdefault(key, replace(
            o, key=key, sequence=sequence, positions=positions,
            breaks=tuple(sorted((n - 1 - j, n - 1 - i, kind) for i, j, kind in o.breaks)),
            reverse_complement=not o.reverse_complement, qualities=o.qualities[::-1]))
    index = _ObservationIndex(list(observations.values()))
    paths, notes = {}, set()
    for seed in sorted(chosen.values(), key=lambda o: (-len(o.sequence), o.sequence, repr(o.positions))):
        key = seed.sequence, tuple(seed.positions)
        if any(_contains(path, key) for path in paths):
            continue
        if len(paths) >= max_paths:
            notes.add("path_limit")
            break
        for sequence, positions, _, _ in _extend_both(
                *key, index, min_overlap, (0, 0, 0), max_paths - len(paths), [], notes):
            paths[sequence, tuple(positions)] = None
            if len(paths) >= max_paths:
                notes.add("path_limit")
                break
    maximal = [p for p in paths if not any(len(q[0]) > len(p[0]) and _contains(q, p) for q in paths)]
    return maximal, notes


def _discover(groups, regions, collector, max_discovery_reads, max_paths, min_overlap):
    chosen, reasons, cache = {}, Counter(), _GeometryCache()
    considered, selected = 0, 0
    for template in groups.templates():
        for _, observations in _observations(template, collector, reasons, regions, cache):
            considered += 1
            if selected >= max_discovery_reads:
                continue
            selected += 1
            for o in observations:
                chosen.setdefault((o.sequence, o.positions), o)
    maximal, notes = _assemble(chosen, max_paths, min_overlap)
    return maximal, dict(target_segments=considered, selected_segments=selected,
                         distinct_observations=len(chosen), max_discovery_reads=max_discovery_reads,
                         max_paths=max_paths, min_overlap=min_overlap,
                         truncated=considered > selected or bool(notes), limitations=sorted(notes),
                         observation_rejections=dict(reasons))


def _support(groups, candidates, regions, collector, header, sample_id, source, *, label="isoform"):
    compatible_key, overlapping_key = "compatible_" + label + "s", "overlapping_" + label + "s"
    ids_key = label + "_ids"
    coverage = {c["id"]: [0] * len(c["sequence"]) for c in candidates}
    supporters, whole, assignments = defaultdict(set), defaultdict(set), Counter()
    whole_cds = defaultdict(set)
    witnesses, reasons, cache = [], Counter(), _GeometryCache()
    lineage = ReadLineage(groups, header)
    labels = CellUmiEvidence(groups, header, sample_id, source, template_loader=groups.template)
    target_fragments, unassigned = 0, 0
    for template in groups.templates():
        observations = list(_observations(template, collector, reasons, regions, cache))
        if not observations:
            continue
        target_fragments += 1
        fragment_matches, segment_rows = None, []
        outside = any(not _relevant(read, regions) for read in template)
        uncertainty = _placement_uncertainty(template)
        outside |= bool(uncertainty)
        for identity, observed in observations:
            compatible, supporting = set(), set()
            for candidate in candidates:
                cid = candidate["id"]
                covered, full = set(), False
                for o in observed:
                    for lo, hi in _matches(candidate["sequence"], candidate["placements"], o):
                        covered.update(range(lo, hi))
                        full |= lo == 0 and hi == len(candidate["sequence"])
                        for i, coding in enumerate(candidate["coding"]):
                            if coding["status"] == "complete_CDS" and lo <= coding["query_start"] and hi >= coding["query_end"]:
                                whole_cds[cid, i].add(identity)
                    if any(_matches(candidate["sequence"], candidate["placements"], o, contained=True)):
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
        assignments[ids, outside] += 1
        unassigned += not ids
        witnesses.append(dict(fragment=list(observations[0][0][:2]), **{compatible_key: list(ids)},
                              competing_placement_outside_targets=outside, placement_uncertainty=uncertainty, segments=segment_rows,
                              records=[dict(id=_record_id(read), contig=read.reference_name,
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
        labels.cache.clear()
        lineage.cache.clear()
    classes = [dict(**{ids_key: list(ids)}, outside_targets=outside, fragments=n)
               for (ids, outside), n in sorted(assignments.items())]
    return dict(scope="complete_declared_input", target_fragments=target_fragments,
                unassigned_target_fragments=unassigned, equivalence_classes=classes,
                observation_rejections=dict(reasons), witnesses=witnesses)



def _placement_uncertainty(template):
    """Known unavailable placements cannot establish unique attribution."""
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
                        reported = _alignment(contig, start, strand == "-", cigar)
                        if not any(r is not read and not r.is_unmapped and _matches_reported(
                                _alignment(r.reference_name, r.reference_start, r.is_reverse, r.cigarstring), reported)
                                   for r in records):
                            reasons.add("missing_SA_placement")
                except ValueError:
                    reasons.add("unparseable_SA_tag")
    return sorted(reasons)
