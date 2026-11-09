"""Bounded original-record acquisition around nominated SV adjacencies."""

from __future__ import annotations

from collections import defaultdict, deque

from .chimeric_alignment import sa_entries
from .default_parameters import SV_BREAKPOINT_WINDOW, SV_MAX_QUERIES, SV_MAX_RECORDS
from .read_identity import sam_record_id, segment_identity
from .rna_path_geometry import interval_contains, merge_intervals

_HOP_WINDOW = 1000  # Nearby SA/mate targets share one indexed query.


def _hop_targets(read):
    """SA and mate locations: retrieval hints only, never observations."""
    targets, notes = [], set()
    if read.has_tag("SA"):
        try:
            targets.extend((contig, position) for contig, position, *_ in sa_entries(read.get_tag("SA")))
        except ValueError:
            notes.add("unparseable_SA_tag")
    if read.is_paired:
        if read.next_reference_id >= 0 and read.next_reference_start >= 0:
            targets.append((read.next_reference_name, read.next_reference_start))
        else:
            notes.add("mate_location_unavailable")
    return targets, notes


def collect_sv_records(bam, regions, references, max_records=SV_MAX_RECORDS, max_queries=SV_MAX_QUERIES,
                       breakpoints=(), breakpoint_window=SV_BREAKPOINT_WINDOW):
    """Fetch both breakpoint neighborhoods and their links before background.

    Parameters
    ----------
    bam : pysam.AlignmentFile
        Indexed alignments.
    regions : iterable of (str, int, int)
        Explicit 0-based half-open search intervals. Exons of every reference
        model are searched too.
    references : iterable of FusionReference
    max_records, max_queries : int
        Resource limits. Reaching one is reported, never silent.
    breakpoints : iterable of FusionBreakpoint
        Each is searched ``breakpoint_window`` bases to either side. Distinct
        windows take turns admitting records, then their observed SA/mate
        links are recovered before other regional/reference-exon queries.

    Returns
    -------
    records : list of pysam.AlignedSegment
        Unmodified records, including competing placements, in a stable order.
    acquisition : dict
        Requested intervals, actual regional/hop queries and limitations.
        Records reachable only through unplaced mates or genome-wide
        realignment are not assessed.
    """
    lengths = dict(zip(bam.references, bam.lengths))
    spans, priority = defaultdict(list), set()
    for breakpoint in breakpoints:
        if breakpoint.contig not in lengths:
            raise ValueError("Breakpoint contig %r is absent from the alignments" % breakpoint.contig)
        start = max(0, breakpoint.position - breakpoint_window)
        end = min(lengths[breakpoint.contig], breakpoint.position + breakpoint_window)
        spans[breakpoint.contig].append((start, end))
        priority.add((breakpoint.contig, start, end, None))
    for contig, start, end in list(regions) + [(r.contig, a, b) for r in references for a, b in r.exons]:
        if (contig not in lengths or type(start) is not int or type(end) is not int
                or not 0 <= start < end <= lengths[contig]):
            raise ValueError("Invalid or unavailable SV search interval: %r" % ((contig, start, end),))
        spans[contig].append((start, end))
    searched = {contig: merge_intervals(intervals) for contig, intervals in spans.items()}
    records, hops, region_queries, requested, limitations = {}, [], [], set(), set()
    covered = defaultdict(list)
    queries = 0

    def collect(pending, interleave=False):
        """Only breakpoint iterators need independent handles and fair turns."""
        targets = defaultdict(set)
        rows = []
        for contig, start, end, fragments in pending:
            row = dict(contig=contig, start=start, end=end, fetched=False, complete=False)
            if fragments is None:
                row['priority'] = interleave
                region_queries.append(row)
            else:
                row['fragments'] = len(fragments)
                hops.append(row)
            rows.append(row)

        def scan(query, row):
            nonlocal queries
            contig, start, end, fragments = query
            if queries >= max_queries or "record_limit" in limitations:
                limitations.add("query_limit" if queries >= max_queries else "record_limit")
                return
            queries += 1
            row['fetched'] = True
            for read in bam.fetch(contig, start, end, multiple_iterators=interleave):
                if read.query_name is None:
                    limitations.add("unnamed_records_skipped")
                    continue
                fragment = segment_identity(read)[:2]
                if fragments is not None and fragment not in fragments:
                    continue
                key = sam_record_id(read)
                if key in records:
                    continue
                if len(records) >= max_records:
                    limitations.add("record_limit")
                    return
                records[key] = read
                hop_targets, notes = _hop_targets(read)
                limitations.update(notes)
                for target_contig, position in hop_targets:
                    if target_contig not in lengths or not 0 <= position < lengths[target_contig]:
                        limitations.add("unavailable_SA_or_mate_contig")
                    elif (target_contig, position, fragment) not in requested:
                        requested.add((target_contig, position, fragment))
                        targets[target_contig, position].add(fragment)
                yield
            row['complete'] = True
            if fragments is None:
                covered[contig] = merge_intervals(covered[contig] + [(start, end)])

        active = deque(scan(query, row) for query, row in zip(pending, rows))
        try:
            while active and "record_limit" not in limitations:
                try:
                    next(active[0])
                except StopIteration:
                    active.popleft()
                else:
                    if interleave:
                        active.rotate(-1)
        finally:
            for iterator in active:
                iterator.close()
        return targets

    def recover(targets):
        while targets:
            pending = []
            for contig, position in sorted(targets):
                # Planned background regions are not evidence of retrieval.
                # Only exhausted, unfiltered queries can make a hop redundant.
                if interval_contains(covered[contig], position):
                    continue
                if pending and pending[-1][0] == contig and position - pending[-1][2] < _HOP_WINDOW:
                    previous = pending[-1]
                    pending[-1] = (contig, previous[1], position + 1, previous[3] | targets[contig, position])
                else:
                    pending.append((contig, position, position + 1, frozenset(targets[contig, position])))
            targets = collect(pending)
            if limitations & {"query_limit", "record_limit"}:
                break

    recover(collect(sorted(priority), interleave=True))
    if not limitations & {"query_limit", "record_limit"}:
        background = [(c, a, b, None) for c in sorted(searched) for a, b in searched[c]
                      if not any(start <= a and b <= end for start, end in covered[c])]
        recover(collect(background))
    return [records[key] for key in sorted(records)], dict(
        searched_regions=[[c, a, b] for c in sorted(searched) for a, b in searched[c]],
        region_queries=region_queries, hop_queries=hops, queries=queries, records=len(records),
        limitations=sorted(limitations),
        scope="indexed_regions_and_observed_SA_or_mate_locations",
        genome_wide_competing_placements_assessed=False)
