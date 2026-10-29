"""Assess every retained SV RNA path against original source evidence.

Sequence matches, exon links and an annotated frame are separate observations.
An exon-linked read matching ordinary shared sequence never transfers its frame
to another RNA branch. Counts are source/read-group scoped, not molecule counts.
"""

from collections import defaultdict
from collections.abc import Iterable, Sequence
from copy import deepcopy
from dataclasses import asdict, dataclass
from hashlib import sha256
import json
from typing import Any, NotRequired, TypedDict

import pysam

from .fusion import FusionReference
from .read_identity import read_group
from .read_collector import ReadCollector
from .read_end_inference import ReadSequenceView
from .transcript_path_frames import TranscriptPathModel, exact_anchor_runs
from .sv_rna_path_support import exact_path_matches, interval_spots, overlap_chain_search, reverse_complement
from .sv_rna_relations import EVENT_CROSSING_RELATIONS


def path_positions(path):
    """Expand observed blocks without filling unplaced RNA from reference."""
    positions = [None] * len(path["sequence"])
    for block in path["blocks"]:
        bases = range(block["reference_start"], block["reference_end"])
        if block["strand"] == "-":
            bases = reversed(bases)
        for q, p in zip(range(block["query_start"], block["query_end"]), bases):
            positions[q] = (block["contig"], p, block["strand"])
    return positions


def coding_intervals(reference):
    """Return CDS genomic blocks in translation order, on either strand."""
    if reference.cds_start is None:
        return []
    offset, intervals = 0, []
    exons = reference.exons if reference.strand == "+" else reference.exons[::-1]
    for start, end in exons:
        left = max(reference.cds_start, offset) - offset
        right = min(reference.cds_end, offset + end - start) - offset
        if left < right:
            intervals.append((start + left, start + right) if reference.strand == "+"
                             else (end - right, end - left))
        offset += end - start
    return intervals


def exon_catalog(references):
    """Index versioned exons and CDS intervals without assuming RNA polarity."""
    exons = defaultdict(lambda: dict(transcripts=set(), coding_intervals=[]))
    starts = []
    for reference in references:
        coding = coding_intervals(reference)
        if coding:
            start = coding[0][0] if reference.strand == "+" else coding[0][1] - 1
            starts.append(dict(transcript_id=reference.transcript_id, contig=reference.contig,
                               strand=reference.strand, position_0based=start))
        for a, b in reference.exons:
            exon = exons[reference.contig, a, b]
            exon["transcripts"].add(reference.transcript_id)
            exon["coding_intervals"].extend((lo, hi, reference.transcript_id)
                                            for lo, hi in coding if lo < b and a < hi)
    return exons, sorted(starts, key=lambda x: x["transcript_id"])


def exon_placements(read, exons, minimum=15):
    """Return actual exon-aligned blocks; skipped introns do not count."""
    if read.is_unmapped:
        return []
    rows = []
    for lo, hi in read.get_blocks():
        for (contig, start, end), exon in exons.items():
            overlap = min(hi, end) - max(lo, start)
            if contig == read.reference_name and overlap >= minimum:
                rows.append(dict(contig=contig, exon_start_0based=start, exon_end_0based=end,
                                 aligned_exon_bases=overlap, transcript_ids=sorted(exon["transcripts"]),
                                 coding_transcript_ids=sorted({tid for a, b, tid in exon["coding_intervals"]
                                                              if min(b, hi) - max(a, lo) >= minimum})))
    return rows


def aligned_exon_clip_contains(record, segment):
    """Require this exon alignment's own query to carry the clipped match."""
    query = record["query_sequence"]
    if not query:
        return False
    start, end = record["query_aligned_interval"]
    for sequence in (segment, reverse_complement(segment)):
        offset = query.find(sequence)
        while offset >= 0:
            eligible = any(left <= offset and offset + len(sequence) <= right
                           for left, right in record.get('usable_query_intervals', [(0, len(query))]))
            if eligible and (offset < start or offset + len(sequence) > end):
                return True
            offset = query.find(sequence, offset + 1)
    return False


RawRead = tuple[str, str, int, str]
Fragment = tuple[str, str | None]
SequenceEvidence = tuple[Fragment, int, str, int | None, int | None]


@dataclass(frozen=True)
class PathEvidenceScope:
    """Inputs that constrain reusable evidence, independent of its counts.

    Exhaustive evidence has no path restriction. Routed evidence is valid only
    for the declared sequences and matching threshold; otherwise a new path or
    shorter match could silently miss original reads excluded by routing.
    """

    references: tuple[FusionReference, ...]
    min_exon_overlap: int
    path_sequences: tuple[str, ...] | None = None
    min_path_match: int | None = None
    sequence_policy: str | None = None

    def validate(self, references: tuple[FusionReference, ...], min_exon_overlap: int,
                 sequences: tuple[str, ...], min_path_match: int, sequence_policy: str | None = None) -> None:
        """Reject reuse that would give an assessment incomplete evidence."""
        if self.sequence_policy is not None and self.sequence_policy != sequence_policy:
            raise ValueError('Path evidence sequence policy differs from assessment')
        if self.references != references or self.min_exon_overlap != min_exon_overlap:
            raise ValueError('Path evidence reference models or exon threshold differ from assessment')
        if self.path_sequences is not None and (
                self.path_sequences != sequences or self.min_path_match != min_path_match):
            raise ValueError('Routed path evidence sequences or match threshold differ from assessment')


class PathEvidence(TypedDict):
    """Original sequences, placements and source-wide accounting for one view."""

    reads: list[SequenceEvidence]
    exon_links: dict[Fragment, list[dict[str, Any]]]
    placements: dict[tuple[Fragment, int], list[dict[str, Any]]]
    records_scanned: int
    fragments_scanned: int
    exon_fragments: int
    coding_exon_fragments: int
    coding_reference_starts: list[dict[str, Any]]
    # Legacy caller-created dictionaries may omit scope. Library collectors
    # always attach it so narrower acquisition is checked before reuse.
    scope: NotRequired[PathEvidenceScope]
    sequence_views: NotRequired[list[dict[str, Any]]]


def collect_path_evidence(bam: pysam.AlignmentFile, references: Sequence[FusionReference],
                          raw_reads: Iterable[RawRead] = (), min_exon_overlap: int = 15,
                          *, read_collector: ReadCollector | None = None) -> PathEvidence:
    """Read a single library, retaining alternative alignments and raw ends.

    ``raw_reads`` supplies (read_group, query_name, mate, original_sequence).
    Mate is 0, 1 or 2. Raw sequences add sequence evidence, not placements.
    The caller must not combine independently processed copies as libraries.
    """
    bam.reset()
    return evidence_from_records(bam.fetch(until_eof=True), references, raw_reads, min_exon_overlap,
                                 read_collector=read_collector)


def _view_sequences(fragment: Fragment, mate: int, view: ReadSequenceView,
                    aligned_start: int | None, aligned_end: int | None) -> Iterable[SequenceEvidence]:
    """Keep separate contiguous chunks under the same original read identity."""
    for start, end in view.usable_intervals:
        yield (fragment, mate, view.original_sequence[start:end].upper(),
               None if aligned_start is None else max(0, min(end - start, aligned_start - start)),
               None if aligned_end is None else max(0, min(end - start, aligned_end - start)))


def _view_record(fragment: Fragment, mate: int, view: ReadSequenceView, alignment: str | None) -> dict[str, Any]:
    """Expose original coordinates and uncertainty without duplicating raw bases."""
    return dict(fragment_id=list(fragment), mate=mate, alignment=alignment,
                original_query_length=len(view.original_sequence), original_sequence=view.original_sequence,
                usable_query_intervals=[list(interval) for interval in view.usable_intervals],
                barrier_intervals=[list(interval) for interval in view.barrier_intervals],
                annotations=[annotation._asdict() for annotation in view.annotations],
                adapter_order_status=view.adapter_order_status,
                upstream_poly_a_length=view.upstream_poly_a_length, upstream_poly_a_status=view.upstream_poly_a_status,
                profile=list(view.profile), reverse=view.reverse, hard_clips=list(view.hard_clips),
                genomic_evidence=None if view.genomic_evidence is None else asdict(view.genomic_evidence))


def evidence_from_records(alignments: Iterable[pysam.AlignedSegment], references: Sequence[FusionReference],
                          raw_reads: Iterable[RawRead] = (), min_exon_overlap: int = 15,
                          *, source_records: int | None = None, source_fragments: int | None = None,
                          read_collector: ReadCollector | None = None) -> PathEvidence:
    """Interpret originals through the same collector for exhaustive and routed use.

    Routed acquisition must supply every potentially matching sequence, all its
    template placements and every alignment to the view's reference exons. Its
    original-library denominators are supplied separately from selected records.
    """
    collector = ReadCollector() if read_collector is None else read_collector
    sequence_views = []
    exons, starts = exon_catalog(references)
    by_contig: dict[str, dict[tuple[str, int, int], dict[str, Any]]] = defaultdict(dict)
    for key, value in exons.items():
        by_contig[key[0]][key] = value
    reads: set[SequenceEvidence] = set()
    records: dict[tuple[Fragment, int], list[dict[str, Any]]] = defaultdict(list)
    exon_links: dict[Fragment, list[dict[str, Any]]] = defaultdict(list)
    fragments, coding, count = set(), set(), 0
    for read in alignments:
        count += 1
        fragment = (read_group(read), read.query_name)
        mate = 1 if read.is_read1 else 2 if read.is_read2 else 0
        fragments.add(fragment)
        view = None
        if read.query_sequence:
            view = (collector.read_sequence_view(read) if collector.sequence_views_enabled else
                    ReadSequenceView(read.query_sequence, None, 0, len(read.query_sequence)))
            reads.update(_view_sequences(fragment, mate, view,
                         None if read.is_unmapped else read.query_alignment_start,
                         None if read.is_unmapped else read.query_alignment_end))
            if view.annotations:
                sequence_views.append(_view_record(fragment, mate, view, sha256(read.to_string().encode()).hexdigest()))
        placements = exon_placements(read, by_contig.get(read.reference_name, {}) if read.reference_name is not None else {},
                                     min_exon_overlap)
        # Keep every original alternative for cited exon links. Record counts
        # never enter the fragment support denominator.
        if placements:
            qualities = read.query_qualities
            lo, hi = read.query_alignment_start, read.query_alignment_end
            row = dict(fragment_id=list(fragment), mate=mate, contig=read.reference_name,
                       reference_start_0based=read.reference_start,
                       reference_end_0based=read.reference_end, cigar=read.cigarstring,
                       mapq=read.mapping_quality, flag=read.flag, reverse=read.is_reverse,
                       secondary=read.is_secondary, supplementary=read.is_supplementary,
                       query_sequence=read.query_sequence, query_qualities=list(qualities) if qualities is not None else None,
                       query_aligned_interval=[lo, hi],
                       left_clip_min_q=min(qualities[:lo]) if qualities is not None and lo else None,
                       right_clip_min_q=min(qualities[hi:]) if qualities is not None and hi < len(qualities) else None,
                       exons=placements)
            if view is not None and view.usable_intervals != ((0, len(view.original_sequence)),):
                row['usable_query_intervals'] = list(view.usable_intervals)
            if row not in exon_links[fragment]:
                exon_links[fragment].append(row)
            if any(p["coding_transcript_ids"] for p in placements):
                coding.add(fragment)
        if read.query_sequence:
            key = (fragment, mate)
            row = dict(fragment_id=list(fragment), mate=mate, contig=read.reference_name,
                       reference_start_0based=read.reference_start, cigar=read.cigartuples,
                       query_sequence=read.query_sequence.upper(), mapq=read.mapping_quality, flag=read.flag,
                       reverse=read.is_reverse, secondary=read.is_secondary,
                       supplementary=read.is_supplementary, unmapped=read.is_unmapped)
            if view is not None and view.usable_intervals != ((0, len(view.original_sequence)),):
                row['usable_query_intervals'] = list(view.usable_intervals)
            if row not in records[key]:
                records[key].append(row)
    for group, name, mate, sequence in raw_reads:
        if type(mate) is not int or mate not in (0, 1, 2) or not isinstance(sequence, str) or not sequence:
            raise ValueError("Raw RNA needs a sequence and a mate number in 0, 1, 2")
        view = (collector.raw_sequence_view(sequence, read_group=group, mate=mate) if collector.sequence_views_enabled else
                ReadSequenceView(sequence, None, 0, len(sequence)))
        reads.update(_view_sequences((group, name), mate, view, None, None))
        if view.annotations:
            sequence_views.append(_view_record((group, name), mate, view, None))
    return dict(reads=sorted(reads, key=repr), exon_links=exon_links, placements=records,
                records_scanned=count if source_records is None else source_records,
                fragments_scanned=len(fragments) if source_fragments is None else source_fragments,
                exon_fragments=len(exon_links), coding_exon_fragments=len(coding),
                coding_reference_starts=starts,
                scope=PathEvidenceScope(tuple(references), min_exon_overlap,
                                        sequence_policy=collector.sequence_policy_fingerprint),
                sequence_views=sequence_views)


def _matched_reads(paths, reads, seed_size):
    """Route originals to paths once, without rescanning a library per path."""
    seeds = defaultdict(set)
    if not paths:
        return defaultdict(list)
    for path in paths:
        sequence = path["sequence"]
        for offset in range(len(sequence) - seed_size + 1):
            seeds[sequence[offset:offset + seed_size]].add(path["path_id"])
    routed = defaultdict(list)
    for read in reads:
        sequence = read[2]
        ids = set()
        for query in (sequence, reverse_complement(sequence)):
            for i in range(len(query) - seed_size + 1):
                ids.update(seeds.get(query[i:i + seed_size], ()))
        for identifier in ids:
            routed[identifier].append(read)
    return routed


def _relevant_sequence_views(paths, evidence, minimum):
    views = evidence.get('sequence_views', [])
    reads = [(tuple(view['fragment_id']), view['mate'], view['original_sequence'].upper(), None, None) for view in views]
    size = min(15, minimum)
    routed = _matched_reads(paths, reads, size)
    fragments = set(evidence['exon_links'])
    for path in paths:
        fragments.update(match.spot for match in exact_path_matches(path['sequence'],
                         routed[path['path_id']], minimum, size))
    return [view for view in views if tuple(view['fragment_id']) in fragments]


def _support(matches, start, end, min_overlap, max_chain_reads):
    if end <= start:
        return dict(interval=[start, end], complete_fragments=[], overlap_chain=[])
    search = overlap_chain_search(matches, start, end, min_overlap=min_overlap, max_reads=max_chain_reads)
    return dict(interval=[start, end], complete_fragments=interval_spots(matches, start, end),
                overlap_chain=search.pop("chain"), chain_search=search)



def _placement_positions(record, start, end, reverse):
    """Map a complete local query interval; clips/insertions remain unplaced."""
    positions = [None] * (end - start)
    q, r = 0, record["reference_start_0based"]
    for operation, length in record["cigar"] or ():
        if operation in (0, 7, 8):
            for offset in range(max(start, q), min(end, q + length)):
                positions[offset - start] = (record["contig"], r + offset - q, "-" if reverse else "+")
            q, r = q + length, r + length
        elif operation in (1, 4):
            q += length
        elif operation in (2, 3):
            r += length
    return positions[::-1] if reverse else positions


def _alternative_placements(window, expected_positions, fragments, placements):
    # These are competing original alignments, not reference transcripts or
    # independent molecules. MAPQ and all provenance are retained for review.
    for fragment in fragments:
        for mate in (0, 1, 2):
            for record in placements.get((fragment, mate), []):
                if record["unmapped"]:
                    continue
                for reverse, segment in ((False, window), (True, reverse_complement(window))):
                    offset = record["query_sequence"].find(segment)
                    while offset >= 0:
                        positions = _placement_positions(record, offset, offset + len(window), reverse)
                        eligible = any(left <= offset and offset + len(window) <= right for left, right in
                                       record.get('usable_query_intervals', [(0, len(record['query_sequence']))]))
                        if all(p is not None for p in positions) and positions != expected_positions:
                            value = dict(sequence=window, positions=positions, original_alignment=record,
                                         query_interval=[offset, offset + len(window)])
                            if not eligible:
                                value['sequence_policy_eligible'] = False
                            identifier = sha256(json.dumps(value, sort_keys=True).encode()).hexdigest()
                            yield dict(placement_id=identifier, **value)
                        offset = record["query_sequence"].find(segment, offset + 1)


def assess_sv_rna_paths(result: dict[str, Any], bam: pysam.AlignmentFile | None,
                        references: Iterable[FusionReference | dict[str, Any]], *,
                        raw_reads: Iterable[RawRead] = (), min_path_match: int = 20,
                        min_exon_overlap: int = 15, junction_flank: int = 12, min_overlap: int = 30,
                        max_chain_reads: int = 12, evidence: PathEvidence | None = None,
                        read_collector: ReadCollector | None = None) -> dict[str, Any]:
    """Assess all retained paths; never select a winner or infer a new frame.

    Returns a JSON-serializable source-specific assessment. Existing annotated
    translations retain their frame evidence. Exploratory ATGs stay separate.
    Exact matching can undercount noisy reads; a missing match is not absence.
    Reused library-collected evidence must have the same reference models and
    exon threshold; routed evidence also pins the path sequences and match
    threshold. Caller-created legacy evidence dictionaries remain accepted.
    """
    collector = ReadCollector() if read_collector is None else read_collector
    for value in (min_path_match, min_exon_overlap, junction_flank, min_overlap, max_chain_reads):
        if type(value) is not int or value < 1:
            raise ValueError("Assessment thresholds must be positive integers")
    reference_models = tuple(FusionReference(**r) if isinstance(r, dict) else r for r in references)
    if any(r.reference_name != result["reference_name"] for r in reference_models):
        raise ValueError("Assessment reference assembly mismatch")
    if evidence is None:
        if bam is None:
            raise ValueError('Assessment requires an original BAM or collected path evidence')
        evidence = collect_path_evidence(bam, reference_models, raw_reads, min_exon_overlap, read_collector=collector)
    paths = result["paths"]
    if 'scope' in evidence:
        evidence['scope'].validate(reference_models, min_exon_overlap,
                                   tuple(path['sequence'] for path in paths), min_path_match, collector.sequence_policy_fingerprint)
    if len({p["path_id"] for p in paths}) != len(paths):
        raise ValueError("Duplicate RNA path identity")
    seed_size = min(15, min_path_match)
    routed = _matched_reads(paths, evidence["reads"], seed_size)
    alignment_ids = {id(record): sha256(json.dumps(record, sort_keys=True).encode()).hexdigest()
                     for records in evidence["exon_links"].values() for record in records}
    record_labels = {digest: "exon:%d" % index for index, digest in enumerate(sorted(set(alignment_ids.values())))}
    alignment_ids = {key: record_labels[digest] for key, digest in alignment_ids.items()}
    models = [TranscriptPathModel(reference) for reference in reference_models]
    rows: list[dict[str, Any]] = []
    alternative_placements: dict[str, dict[str, Any]] = {}
    exon_alignments: dict[str, dict[str, Any]] = {}
    for path in paths:
        sequence = path["sequence"]
        matches = exact_path_matches(sequence, routed[path["path_id"]], min_path_match, seed_size)
        positions = path_positions(path)
        anchors = [dict(transcript_id=model.reference.transcript_id,
                        annotation=model.reference.annotation, query_interval=[a, b],
                        transcript_interval=[offset, offset + b - a],
                        coding_bases=max(0, min(offset + b - a, model.reference.cds_end or 0)
                                         - max(offset, model.reference.cds_start or 0)))
                   for model in models for a, b, offset in exact_anchor_runs(sequence, positions, model)
                   if b - a >= min_path_match]
        joins = []
        for index, junction in enumerate(path["junctions"]):
            left, right = junction["query_interval"]
            start, end = max(0, left - junction_flank), min(len(sequence), right + junction_flank)
            support = _support(matches, start, end, min_overlap, max_chain_reads)
            for alternative in _alternative_placements(
                    sequence[start:end], positions[start:end], support["complete_fragments"], evidence["placements"]):
                key = alternative["placement_id"]
                if key not in alternative_placements:
                    alternative_placements[key] = alternative
            joins.append(dict(junction_index=index, relation=junction["relation"],
                              query_interval=[left, right], window_sequence=sequence[start:end],
                              direct_support=deepcopy(junction["direct_support"]),
                              direct_observations=deepcopy(junction.get("direct_observations", [])),
                              exact_sequence_support=support,
                              exon_linked_fragments=[fragment for fragment in support["complete_fragments"]
                                                     if fragment in evidence["exon_links"]]))
        links: dict[tuple[Fragment, int, str], dict[str, Any]] = {}
        for match in matches:
            for record in evidence["exon_links"].get(match.spot, []):
                same_end = match.mate == record["mate"]
                record_id = alignment_ids[id(record)]
                exon_alignments.setdefault(record_id, record)
                key = match.spot, match.mate, record_id
                link = links.setdefault(key, dict(fragment_id=list(match.spot), mate=match.mate,
                    same_end_exon_alignment=same_end, opposite_mate_exon_alignment={match.mate, record["mate"]} == {1, 2},
                    exon_alignment_id=record_id, matches=[]))
                link["matches"].append(dict(path_interval=[match.start, match.end],
                    same_end_exon_clip_contains_match=same_end and aligned_exon_clip_contains(
                        record, sequence[match.start:match.end])))
        proteins = []
        for candidate in path.get("translations", []):
            linked = bool(set(candidate.get("departure_relations", [])) & set(EVENT_CROSSING_RELATIONS))
            # Supplied-fusion translations already require a coding junction.
            linked |= result.get("schema", "").startswith("isovar.fusion_rna.")
            proteins.append(dict(candidate=deepcopy(candidate), event_linked=linked,
                                 completeness="complete" if candidate.get("complete_5prime") and
                                 candidate["ends_with_stop_codon"] else "partial",
                                 sequence_support=_support(matches, candidate["translation_start"],
                                                           candidate["translation_end"], min_overlap, max_chain_reads)))
        exploratory = []
        for candidate in path.get("exploratory_orfs", {}).get("candidates", []):
            a, b = candidate["query_interval"]
            exploratory.append(dict(candidate=deepcopy(candidate),
                                    sequence_support=_support(matches, a, b, min_overlap, max_chain_reads),
                                    initiation_observed=False, translation_observed=False))
        rows.append(dict(path_id=path["path_id"], sequence=sequence, blocks=deepcopy(path["blocks"]),
                         junctions=joins, exon_links=list(links.values()), exact_exon_anchors=anchors,
                         frame_status=path["frame_status"], proteins=proteins, exploratory_orfs=exploratory,
                         local_branch_support=deepcopy(path.get("event_local_support")),
                         full_sequence_support=_support(matches, 0, len(sequence), min_overlap, max_chain_reads),
                         exact_matches=[dict(fragment_id=list(m.spot), mate=m.mate,
                                             path_interval=[m.start, m.end], clipped=m.clipped) for m in matches]))
    return dict(schema="isovar.sv_rna_assessment.v1", event_id=result["event_id"],
                reference_name=result["reference_name"], sample_id=result["sample_id"],
                source=result.get("source", result.get("provenance", {}).get("source")),
                event_provenance=deepcopy(result.get("event_provenance", result.get("provenance"))),
                paths=rows, records_scanned=evidence["records_scanned"],
                fragments_scanned=evidence["fragments_scanned"],
                exon_fragments=evidence["exon_fragments"], coding_exon_fragments=evidence["coding_exon_fragments"],
                coding_reference_starts=evidence["coding_reference_starts"],
                reference_models=[asdict(r) for r in reference_models],
                exon_alignments=exon_alignments, read_sequence_views=_relevant_sequence_views(paths, evidence, min_path_match),
                alternative_placements=sorted(alternative_placements.values(), key=lambda r: r["placement_id"]),
                competing_annotated_junctions=deepcopy(result.get("competing_annotated_junctions", [])),
                limitations=deepcopy(result.get("limitations", result.get("reasons", []))),
                parameters=dict(min_path_match=min_path_match, min_exon_overlap=min_exon_overlap,
                                junction_flank=junction_flank, min_overlap=min_overlap,
                                max_chain_reads=max_chain_reads),
                sequence_policy_fingerprint=collector.sequence_policy_fingerprint,
                interpretation="sequence_and_frame_evidence_not_observed_translation",
                translation_observed=False, somatic_causation_proven=False)
